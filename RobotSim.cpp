#include <iostream>
#include <fstream>
#include <windows.h>
#include <chrono>
#include <thread>

// zgodne z danymi katalogowymi IRB 930 1.05_0.3

const double globalZoffset = -0.3;

const double arm1 = 0.56;
const double arm2 = 0.49;
const double arm3 = 0.855 - 0.275;

const double fUpLimit = 145.0 / 180.0 * 3.14159265359;
const double fLowLimit = 0.01 / 180.0 * 3.14159265359;
const double hUpLimit = 0.275 + globalZoffset;
const double hLowLimit = hUpLimit - arm3;

const double minArmDist = (arm1 + arm2 * cos(fUpLimit)) * (arm1 + arm2 * cos(fUpLimit)) + arm2 * sin(fUpLimit) * arm2 * sin(fUpLimit);
const double maxArmDist = (arm1 + arm2 * cos(fLowLimit)) * (arm1 + arm2 * cos(fLowLimit)) + arm2 * sin(fLowLimit) * arm2 * sin(fLowLimit);

const double SumSqrArms = arm1 * arm1 + arm2 * arm2;
const double DifSqrArms = arm1 * arm1 - arm2 * arm2;
const double SqrSumArms = (arm1 + arm2) * (arm1 + arm2);
const double SqrDifArms = (arm1 - arm2) * (arm1 - arm2);
const double mapConst1 = 2 * arm1 * arm1 + 2 * arm2 * arm2;
const double mapConst2 = (arm1 * arm1 - arm2 * arm2) * (arm1 * arm1 - arm2 * arm2);

double massCoeffs[6];
double ELSol[9];

struct projectivePlane {
	// 1024 x 683
	double cdist;	// camera distance from global origin
	double cangle[2];	// two spherical angles
	double normal[3];
	double origin[3];	// local origin
	double camera[3];	// camera in global coordinates
	double jacobian[6]; // { du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz }

};

struct joint {

	double pos[3];		// joint position in global basis
	double offset[3];	// joint offset in local basis

	// drawing data
	double vert[8][3];		// joint vertex offsets from pos
	double screen[8][2];	// joint vertex offsets mapped onto screen

	// state variables
	double angle;	// joint angle
	double vel;		// joint velocity
	double drive;	// driving force
	double force;	// total active force
	double reaction;	// joint motion change resistance
	double acc;
	double mass; // joint motion change resistance at rest
	double pressure; // frition normal force

	// control variables
	double PIDsum;	// integral part of PID force regulation
	double prevdrive; // used for ramp damping
	double kp;	// 0 - vel, 1 - acc | kp
	double ki;	// 0 - vel, 1 - acc | kp * Ts / Ti

	double basis[9];	// forward transfrom from global to local

};

void joint_calc(joint* c, joint* n) {

	n->pos[0] = c->offset[0] * c->basis[0] + c->offset[1] * c->basis[3] + c->offset[2] * c->basis[6] + c->pos[0];
	n->pos[1] = c->offset[0] * c->basis[1] + c->offset[1] * c->basis[4] + c->offset[2] * c->basis[7] + c->pos[1];
	n->pos[2] = c->offset[0] * c->basis[2] + c->offset[1] * c->basis[5] + c->offset[2] * c->basis[8] + c->pos[2];
	double C = cos(c->angle);
	double S = sin(c->angle);

	// 5: x=Cx'+Sy'	y=-Sx'+Cy'	z=-z'
	n->basis[0] = c->basis[0] * C + c->basis[3] * S;
	n->basis[1] = c->basis[1] * C + c->basis[4] * S;
	n->basis[2] = c->basis[2] * C + c->basis[5] * S;

	n->basis[3] = -c->basis[0] * S + c->basis[3] * C;
	n->basis[4] = -c->basis[1] * S + c->basis[4] * C;
	n->basis[5] = -c->basis[2] * S + c->basis[5] * C;

	n->basis[6] = c->basis[6];
	n->basis[7] = c->basis[7];
	n->basis[8] = c->basis[8];

}

void projectivePlane_update(projectivePlane* p) {

	// temporarily, camera only looks at 0 
	// if updated, then coordinate functions change

	p->normal[0] = cos(p->cangle[0]) * cos(p->cangle[1]);
	p->normal[1] = sin(p->cangle[0]) * cos(p->cangle[1]);
	p->normal[2] = sin(p->cangle[1]);

	p->jacobian[0] = -sin(p->cangle[0]);
	p->jacobian[1] = cos(p->cangle[0]);
	p->jacobian[2] = 0;
	p->jacobian[3] = cos(p->cangle[0]) * sin(p->cangle[1]);
	p->jacobian[4] = sin(p->cangle[0]) * sin(p->cangle[1]);
	p->jacobian[5] = -cos(p->cangle[1]);

	// origin is 2 global units towards global origin starting from camera,
	// then 1.5 global units in direction of - du,
	// then 1 global units in the direction of -dv
	// assuming 3x2 global unit plane, camera is now in center
	p->camera[0] = p->cdist * p->normal[0];
	p->camera[1] = p->cdist * p->normal[1];
	p->camera[2] = p->cdist * p->normal[2];

	p->origin[0] = p->camera[0] - p->normal[0] - p->normal[0];
	p->origin[1] = p->camera[1] - p->normal[1] - p->normal[1];
	p->origin[2] = p->camera[2] - p->normal[2] - p->normal[2];

	p->origin[0] -= 1.5 * p->jacobian[0];
	p->origin[1] -= 1.5 * p->jacobian[1];

	p->origin[0] -= p->jacobian[3];
	p->origin[1] -= p->jacobian[4];
	p->origin[2] -= p->jacobian[5];

	// scaling jacobian to get basis in terms of pixels on the screen
	p->jacobian[0] *= 341.0;
	p->jacobian[1] *= 341.0;
	p->jacobian[3] *= 341.0;
	p->jacobian[4] *= 341.0;
	p->jacobian[5] *= 341.0;
}

bool projectivePlane_move(projectivePlane* p) {
	bool flag = 0;
	if (GetKeyState(VK_LEFT) & 0x8000) {
		p->cangle[0] += 0.03;
		flag = 1;
	}
	if (GetKeyState(VK_RIGHT) & 0x8000) {
		p->cangle[0] -= 0.03;
		flag = 1;
	}
	if (GetKeyState(VK_UP) & 0x8000) {
		p->cangle[1] += 0.03;
		flag = 1;
	}
	if (GetKeyState(VK_DOWN) & 0x8000) {
		p->cangle[1] -= 0.03;
		flag = 1;
	}
	if (GetKeyState('W') & 0x8000) {
		if (p->cdist > 0.05) {
			p->cdist -= 0.05;
			flag = 1;
		}
	}
	if (GetKeyState('S') & 0x8000) {
		p->cdist += 0.05;
		flag = 1;
	}
	if (GetKeyState('A') & 0x8000) {
		p->cangle[0] = 0;
		p->cangle[1] = 0;
		p->cdist = 3;
		flag = 1;
	}
	return flag;
}

void projectivePlane_map_fromR3(double* x, projectivePlane* pp, double* y) {
	double l[3];
	l[0] = pp->camera[0] - x[0];
	l[1] = pp->camera[1] - x[1];
	l[2] = pp->camera[2] - x[2];
	double d = (pp->origin[0] - x[0]) * pp->normal[0];
	d += (pp->origin[1] - x[1]) * pp->normal[1];
	d += (pp->origin[2] - x[2]) * pp->normal[2];
	d = d / (l[0] * pp->normal[0] + l[1] * pp->normal[1] + l[2] * pp->normal[2]);
	d = abs(d);
	l[0] *= d;
	l[1] *= d;
	l[2] *= d;
	l[0] += x[0];
	l[1] += x[1];
	l[2] += x[2];
	// l is a point on the projection plane
	// identify it with a vector on the plane
	l[0] = l[0] - pp->origin[0];
	l[1] = l[1] - pp->origin[1];
	l[2] = l[2] - pp->origin[2];
	// du = du/dx dx + du/dy dy + du/dz dz
	y[0] = l[0] * pp->jacobian[0] + l[1] * pp->jacobian[1];
	y[1] = l[0] * pp->jacobian[3] + l[1] * pp->jacobian[4] + l[2] * pp->jacobian[5];
	return;
}

void screen_drawLine(double* a, double* b, int* frame, int color) {
	int x = a[0];
	int y = a[1];
	x = x << 18;	// defines start
	y = y << 18;
	double temp1 = (b[0] - a[0]);
	double temp2 = (b[1] - a[1]);
	double temp = temp1 * temp1 + temp2 * temp2;
	temp = sqrt(temp);
	temp1 = 262144.0 / temp;	// shifts velocity to correct constant point
	int v1 = (int)((b[0] - a[0]) * temp1);
	int v2 = (int)(temp2 * temp1);

	double i = 0;
	while ((i < temp)) {

		// 2^18 * 1024		2^18 * 683		check if in bounds
		if ((x < 0x10000000) && (y < 179044352) && (x > 0) && (y > 0)) {
			frame[(x >> 18) + ((y >> 8) & 0xFFFFFC00)] = color; // y must be an integer - erase fractional part
			frame[(x >> 18) + ((y >> 8) & 0xFFFFFC00) + 1024] = color;
		}
		else {
			x += v1 << 3;
			y += v2 << 3;
			i += 8;
		}
		x += v1;
		y += v2;
		i++;
	}
}

void screen_drawJoint(joint* j, int* frame, projectivePlane* p, int color) {

	double temp[3];
	for (int i = 0; i < 8; i++) {
		temp[0] = j->vert[i][0] * j->basis[0] + j->vert[i][1] * j->basis[3] + j->vert[i][2] * j->basis[6] + j->pos[0];
		temp[1] = j->vert[i][0] * j->basis[1] + j->vert[i][1] * j->basis[4] + j->vert[i][2] * j->basis[7] + j->pos[1];
		temp[2] = j->vert[i][0] * j->basis[2] + j->vert[i][1] * j->basis[5] + j->vert[i][2] * j->basis[8] + j->pos[2];
		projectivePlane_map_fromR3(temp, p, j->screen[i]);
	}

	screen_drawLine(j->screen[0], j->screen[4], frame, color);
	screen_drawLine(j->screen[1], j->screen[5], frame, color);
	screen_drawLine(j->screen[2], j->screen[6], frame, color);
	screen_drawLine(j->screen[3], j->screen[7], frame, color);

	screen_drawLine(j->screen[0], j->screen[2], frame, color);
	screen_drawLine(j->screen[0], j->screen[3], frame, color);
	screen_drawLine(j->screen[1], j->screen[2], frame, color);
	screen_drawLine(j->screen[1], j->screen[3], frame, color);

	screen_drawLine(j->screen[4], j->screen[6], frame, color);
	screen_drawLine(j->screen[4], j->screen[7], frame, color);
	screen_drawLine(j->screen[5], j->screen[6], frame, color);
	screen_drawLine(j->screen[5], j->screen[7], frame, color);

}

double friction(double vel, double fdrive, double fpress, double cStatic, double cCoulomb, double cViscous, double InvSqrRefVel) {

	double sign = 1.0;
	fpress = abs(fpress);

	if (vel == 0) {
		if (fdrive < 0) {
			fdrive = -fdrive;
			sign = -1.0;
		}
		if (cStatic * fpress > fdrive) return sign * fdrive;
		return sign * cStatic * fpress;
	}

	if (vel < 0) {
		vel = -vel;
		sign = -1.0;
	}

	return sign * fpress * (cCoulomb + (cStatic - cCoulomb) * exp(-vel * vel * InvSqrRefVel) + cViscous * vel);

}

void initRobot(joint* base, joint* z1, joint* z2, joint* z3, joint* load) {

	memset(base, 0, sizeof(joint));
	memset(z1, 0, sizeof(joint));
	memset(z2, 0, sizeof(joint));
	memset(z3, 0, sizeof(joint));
	memset(load, 0, sizeof(joint));

	base->pos[2] = globalZoffset;
	base->basis[0] = 1;
	base->basis[4] = 1;
	base->basis[8] = 1;

	base->angle = -0.2;

	base->offset[2] = 0.275;
	base->vert[0][0] = .11;
	base->vert[0][1] = -.11;
	base->vert[1][0] = -.11;
	base->vert[1][1] = .11;
	base->vert[3][0] = .11;
	base->vert[3][1] = .11;
	base->vert[2][0] = -.11;
	base->vert[2][1] = -.11;
	for (int i = 0; i < 4; i++) {
		base->vert[i + 4][0] = base->vert[i][0];
		base->vert[i + 4][1] = base->vert[i][1];
		base->vert[i + 4][2] = base->vert[i][2] + 0.275;
	}

	z1->angle = 0.6;
	z1->offset[0] = arm1;
	z1->offset[2] = 0.0;
	z1->vert[0][0] = arm1;
	z1->vert[0][1] = -.1;
	z1->vert[1][0] = -.1;
	z1->vert[1][1] = .1;
	z1->vert[3][0] = arm1;
	z1->vert[3][1] = .1;
	z1->vert[2][0] = -.1;
	z1->vert[2][1] = -.1;
	for (int i = 0; i < 4; i++) {
		z1->vert[i + 4][0] = z1->vert[i][0];
		z1->vert[i + 4][1] = z1->vert[i][1];
		z1->vert[i + 4][2] = z1->vert[i][2] + 0.1;
	}

	z2->offset[0] = arm2;
	z2->offset[2] = -0.2;
	z2->vert[0][0] = arm2 + 0.085;
	z2->vert[0][1] = -.1;
	z2->vert[1][0] = -.1;
	z2->vert[1][1] = .1;
	z2->vert[3][0] = arm2 + 0.085;
	z2->vert[3][1] = .1;
	z2->vert[2][0] = -.1;
	z2->vert[2][1] = -.1;
	for (int i = 0; i < 4; i++) {
		z2->vert[i + 4][0] = z2->vert[i][0];
		z2->vert[i + 4][1] = z2->vert[i][1];
		z2->vert[i + 4][2] = z2->vert[i][2] + 0.3;
	}
	z2->vert[4][2] -= .1;
	z2->vert[7][2] -= .1;
	z2->vert[1][2] -= .1;
	z2->vert[2][2] -= .1;

	z3->vert[0][0] = .02;
	z3->vert[0][1] = -.02;
	z3->vert[1][0] = -.02;
	z3->vert[1][1] = .02;
	z3->vert[3][0] = .02;
	z3->vert[3][1] = .02;
	z3->vert[2][0] = -.02;
	z3->vert[2][1] = -.02;
	for (int i = 0; i < 4; i++) {
		z3->vert[i + 4][0] = z3->vert[i][0];
		z3->vert[i + 4][1] = z3->vert[i][1];
		z3->vert[i + 4][2] = z3->vert[i][2] + arm3;
	}

	z3->offset[2] = -0.3;
	load->vert[0][0] = .15;
	load->vert[0][1] = -.15;
	load->vert[1][0] = -.15;
	load->vert[1][1] = .15;
	load->vert[3][0] = .15;
	load->vert[3][1] = .15;
	load->vert[2][0] = -.15;
	load->vert[2][1] = -.15;
	for (int i = 0; i < 4; i++) {
		load->vert[i + 4][0] = load->vert[i][0];
		load->vert[i + 4][1] = load->vert[i][1];
		load->vert[i + 4][2] = load->vert[i][2] + 0.3;
	}

	joint_calc(base, z1);
	joint_calc(z1, z2);
	joint_calc(z2, z3);
	joint_calc(z3, load);

	load->mass = 10.0;
	base->mass = 15.0;
	z3->mass = 11.0;
	z2->mass = 25.0;
	z1->mass = 15.0;
	z3->mass += load->mass;

	z2->PIDsum = z2->mass * 9.81; // standby with equal force

	base->pressure = 9.81 * (load->mass + z2->mass + z1->mass);
	z1->pressure = 9.81 * (load->mass + z2->mass);
	z2->pressure = 10.0;

	massCoeffs[0] = arm1 * arm1 / 3.0 * base->mass + SumSqrArms * z2->mass + (arm2 * arm2 / 3.0 + arm1 * arm1) * z1->mass;
	massCoeffs[1] = 2.0 * arm1 * arm2 * z2->mass + arm1 * arm2 * z1->mass;
	massCoeffs[2] = arm2 * arm2 * z2->mass + arm2 * arm2 / 3.0 * z1->mass;
	massCoeffs[3] = arm1 * arm2 * z2->mass + arm1 * arm2 * 0.5 * z1->mass;
	massCoeffs[4] = arm2 * arm2 * z2->mass + arm2 * arm2 / 3.0 * z1->mass;
	massCoeffs[5] = z2->mass * 0.5;

	base->kp = 2000.0;
	base->ki = base->kp / 600.0 / 2.0;

	z1->kp = 2000.0;
	z1->ki = z1->kp / 600.0 / 2.0;

	z2->kp = 1000.0;
	z2->ki = z2->kp / 600.0 / 2.0;

}

int main() {

	HDC hdc = GetDC(GetConsoleWindow());
	HDC buf = CreateCompatibleDC(hdc);

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	projectivePlane p;
	p.cangle[0] = 4.5;
	p.cangle[1] = 1.1;
	p.cdist = 3;

	joint base;
	joint z1;
	joint z2;
	joint z3;
	joint load;

	initRobot(&base, &z1, &z2, &z3, &load);

	bool flagDraw = 0;
	bool flagPlaneUpdt = 1;
	bool flagCalculate = 0;

	double dist = 0.0;// , disttrav = 0.0, disttabs = 0.0;

	double r = 0.0;		// jacobians for velocity transform
	double drdx, drdy, dadx, dady;
	double dfdr, dtdr; // dpdt = 0, dfdt = 1

	double v = 0.0, vi = 0;	// x cartesian velocities
	double u = 0.0, ui = 0;	// y
	double w = 0.0, wi = 0;	// z

	double SetVel = 1.0;

	double distKp = 20.0;
	double distKi = distKp / 600.0 / 1.0;

	double vel_t = 0.0; // -> base.angle
	double vel_f = 0.0; // -> z1.angle

	double px = 0.0;	// cartesian destination point
	double py = 0.0;
	double pz = 0.0;

	dist = (px - z3.pos[0]) * (px - z3.pos[0]) + (py - z3.pos[1]) * (py - z3.pos[1])
		+ (pz - z3.pos[2]) * (pz - z3.pos[2]);
	dist = sqrt(dist);

	bool antiwindupz1 = 1;
	bool antiwindupz2 = 1;

	int* frameBuffer = (int*)calloc(1024 * 683, 4);
	HBITMAP hbitmap;

	std::ofstream txtdata;
	double processdata[3];
	memset(processdata, 0, sizeof(processdata));
	txtdata.open("last_execution_data.txt");

	double timeStamp = 0;

	txtdata << "t\t" << "X\t" << "Y\t" << "Z\t" << "F_t\t" << "F_f\t" << "F_h\t" << "V_t\t" << "V_f\t" << "V_h\n";

	char exit = 0;
	bool breaking = 0;

	while (1) {//flagCalculate) {
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		if (flagCalculate) {
			flagDraw = 1;

			for (int i = 0; i < 10; i++) {

				joint_calc(&base, &z1);
				joint_calc(&z1, &z2);
				joint_calc(&z2, &z3);
				joint_calc(&z3, &load);

				// velocity regulation

				dist = (px - z3.pos[0]) * (px - z3.pos[0]) + (py - z3.pos[1]) * (py - z3.pos[1])
					+ (pz - z3.pos[2]) * (pz - z3.pos[2]);
				dist = sqrt(dist);

				// switching to position regulation near desination

				if (dist < SetVel / distKp) breaking = 1;

				if (breaking) {

					v = distKp * (px - z3.pos[0]) + vi;
					vi += distKi * (px - z3.pos[0]);

					u = distKp * (py - z3.pos[1]) + ui;
					ui += distKi * (py - z3.pos[1]);

					w = distKp * (pz - z3.pos[2]) + wi;
					wi += distKi * (pz - z3.pos[2]);

					if (dist < 0.001) {
						i = 11;
						flagCalculate = 0;
						breaking = 0;
					}
				}
				else {

					v = SetVel * (px - z3.pos[0]) / dist;
					u = SetVel * (py - z3.pos[1]) / dist;
					w = SetVel * (pz - z3.pos[2]) / dist;
				}

				// transform cartesian velocity into joint velocity

				r = load.pos[0] * load.pos[0] + load.pos[1] * load.pos[1];

				dadx = -z3.pos[1] / r;
				dady = z3.pos[0] / r;

				dtdr = (r - SqrSumArms) * (SqrDifArms - r);

				dfdr = dtdr;
				dfdr = sqrt(dfdr);

				dtdr = r * dfdr;
				dtdr = (r - DifSqrArms) / dtdr;

				r = sqrt(r);

				dtdr *= r;
				dfdr = -2.0 * r / dfdr;

				drdx = dady * r;
				drdy = -dadx * r;

				vel_t = v * (drdx * dtdr + dadx) + u * (drdy * dtdr + dady);
				vel_f = v * drdx * dfdr + u * drdy * dfdr;

				// displacment

				base.angle += 0.00166666 * base.vel;

				z1.angle += 0.00166666 * z1.vel;
				if (z1.angle > fUpLimit or z1.angle < fLowLimit) {
					z1.angle -= 0.00166666 * z1.vel;
					antiwindupz1 = 0;
				}

				z2.offset[2] += 0.00166666 * z2.vel;
				if (z3.pos[2] < hLowLimit or z3.pos[2] > hUpLimit) {
					z2.offset[2] -= 0.00166666 * z2.vel;
					antiwindupz2 = 0;
				}

				// PI Regulation with antiwindup and ramp damping

				base.drive = base.kp * (vel_t - base.vel) + base.PIDsum;
				base.PIDsum += base.ki * (vel_t - base.vel);
				if (base.drive - base.prevdrive > 50.0) base.drive = base.prevdrive + 50.0;
				if (base.drive - base.prevdrive < -50.0) base.drive = base.prevdrive - 50.0;

				z1.drive = z1.kp * (vel_f - z1.vel) + z1.PIDsum;
				if (antiwindupz1) z1.PIDsum += z1.ki * (vel_f - z1.vel);
				if (z1.drive - z1.prevdrive > 50.0) z1.drive = z1.prevdrive + 50.0;
				if (z1.drive - z1.prevdrive < -50.0) z1.drive = z1.prevdrive - 50.0;

				z2.drive = z2.kp * (w - z2.vel) + z2.PIDsum;
				if (antiwindupz2) z2.PIDsum += z2.ki * (w - z2.vel);
				if (z2.drive - z2.prevdrive > 50.0) z2.drive = z2.prevdrive + 50.0;
				if (z2.drive - z2.prevdrive < -50.0) z2.drive = z2.prevdrive - 50.0;

				antiwindupz1 = 1;
				antiwindupz2 = 1;

				base.prevdrive = base.drive;
				z1.prevdrive = z1.drive;
				z2.prevdrive = z2.drive;

				// preparing equations

				ELSol[5] = massCoeffs[0] + massCoeffs[1] * cos(z1.angle);
				ELSol[6] = massCoeffs[2] + massCoeffs[3] * cos(z1.angle);
				ELSol[7] = ELSol[6];
				ELSol[8] = massCoeffs[4];

				ELSol[0] = ELSol[5] * ELSol[8] - ELSol[6] * ELSol[7];
				ELSol[0] = 1.0 / ELSol[0];

				ELSol[1] = ELSol[8] * ELSol[0];
				ELSol[2] = -ELSol[7] * ELSol[0];
				ELSol[3] = ELSol[2];
				ELSol[4] = ELSol[5] * ELSol[0];

				base.force = base.drive + (massCoeffs[3] * z1.vel + massCoeffs[1] * base.vel) * z1.vel * sin(z1.angle);
				z1.force = z1.drive - 0.5 * massCoeffs[1] * base.vel * base.vel * sin(z1.angle);
				z2.force = z2.drive - 9.81 * z3.mass;

				// friction handling - friction is a function of force only when V=0

				if (z1.vel == 0) {
					z1.reaction = z1.force - ELSol[7] * base.acc + ELSol[8] * z1.acc;
				}
				if (z2.vel == 0) {
					z2.reaction = z2.force - z2.mass * z2.acc;
				}

				z2.force -= friction(z2.vel, z2.reaction, z2.pressure, 0.3, 0.2, 0.02, 10000.0);

				// isolating drive - friction for preassure forces
				z2.reaction = z2.force + 9.81 * z3.mass;

				// isolating drive - friction to add to base force
				z1.reaction = friction(z1.vel, z1.reaction, z1.pressure + z2.reaction, 0.25, 0.2, 0.02, 10000.0);

				base.force = base.force + z1.drive - z1.reaction;
				if (base.vel == 0) {
					base.reaction = base.force - ELSol[5] * base.acc + ELSol[6] * z1.acc;
				}

				// adding friction
				z1.force -= z1.reaction;
				base.force -= friction(base.vel, base.reaction, base.pressure + z2.reaction, 0.25, 0.2, 0.02, 10000.0);

				base.acc = base.force * ELSol[4] - z1.force * ELSol[2];
				z1.acc = -base.force * ELSol[3] + z1.force * ELSol[1];
				z2.acc = z2.force / z2.mass;

				// calculating velocity

				base.vel = base.vel + 0.00166666 * base.acc;
				z1.vel = z1.vel + 0.00166666 * z1.acc;
				z2.vel = z2.vel + 0.00166666 * z2.acc;

				// printing to file

				if (i == 0) {

					txtdata << timeStamp << "\t";
					txtdata << z3.pos[0] << "\t" << z3.pos[1] << "\t" << z3.pos[2] - globalZoffset << "\t";
					txtdata << base.drive << "\t" << z1.drive << "\t" << z2.drive << "\t";
					txtdata << base.vel << "\t" << z1.vel << "\t" << z2.vel << "\n";
				}

				timeStamp += 1.0 / 600.0;
			}
		}

		if (flagPlaneUpdt) {
			flagPlaneUpdt = 0;
			flagDraw = 1;

			projectivePlane_update(&p);
		}

		if (flagDraw) {
			flagDraw = 0;

			screen_drawJoint(&base, frameBuffer, &p, 0xFFFFFF);
			screen_drawJoint(&z1, frameBuffer, &p, 0xFFFF00);
			screen_drawJoint(&z2, frameBuffer, &p, 0xFF00FF);
			screen_drawJoint(&z3, frameBuffer, &p, 0x00FFFF);
			if (load.mass != 0) screen_drawJoint(&load, frameBuffer, &p, 0xA0A0A0);

			hbitmap = CreateBitmap(1024, 683, 1, 32, (void*)frameBuffer);
			SelectObject(buf, hbitmap);
			BitBlt(hdc, 0, 0, 1024, 683, buf, 0, 0, SRCCOPY);

			DeleteObject(hbitmap);

			memset(frameBuffer, 0, 4 * 1024 * 683);
		}

		if (not flagCalculate) {

			std::cout << "Exit? [y/n]\n";
			std::cin >> exit;

			if (exit == 121) { // "y"
				free(frameBuffer);
				txtdata.close();

				return 0;
			}

			std::cout << "Current z3 Position\nX = " << z3.pos[0] << "\tY = " << z3.pos[1] << "\tZ = " << z3.pos[2] - globalZoffset;
			std::cout << "\nX :=\t";
			std::cin >> px;
			std::cout << "Y :=\t";
			std::cin >> py;
			std::cout << "Z :=\t";
			std::cin >> pz;
			pz += globalZoffset;
			std::cout << "Velocity [m/s] :=\t";
			std::cin >> SetVel;

			// zerowanie zmiennych pętli obliczeniowej

			vi = 0;
			ui = 0;
			wi = 0;
			dist = (px - z3.pos[0]) * (px - z3.pos[0]) + (py - z3.pos[1]) * (py - z3.pos[1])
				+ (pz - z3.pos[2]) * (pz - z3.pos[2]);
			dist = sqrt(dist);

			if (px * px + py * py < maxArmDist and px * px + py * py > minArmDist and pz < hUpLimit and pz > hLowLimit) {
				flagCalculate = 1;
				system("CLS");
			}

			else {
				std::cout << "\nOut of Bounds\n";
			}
		}

		flagPlaneUpdt = projectivePlane_move(&p);

		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		std::this_thread::sleep_for(std::chrono::microseconds(8333) - (t2 - t1));
	}

}