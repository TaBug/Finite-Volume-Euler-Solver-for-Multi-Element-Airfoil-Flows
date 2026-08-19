#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include "tools.h"
using namespace std;

struct structFlux {
	vector<double> F;
	double s_mag;
};

/* Each flux exists once, as a "Core" routine that reads four doubles per state
 * and writes four doubles of flux through a caller-supplied pointer, returning
 * the wave speed. Nothing here allocates.
 *
 * The vector-taking overloads below are thin wrappers over those cores, kept so
 * getRes.h and test_flux.cpp compile unchanged. They allocate a structFlux per
 * call, which is why the second-order driver calls the cores directly - at
 * 3313 faces x 2 passes x 2 RK stages that allocation was the single largest
 * cost in the solve. Because there is only one implementation, the two entry
 * points cannot drift apart.
 *
 * F must have room for 4 doubles. n is the unit normal, passed as two scalars.
 */

inline double roeCore(const double* UL, const double* UR, double gamma, double nx, double ny, double* F) {
	const double gmi = gamma - 1.0;

	// process left state variables
	double rhoL = UL[0];                        // density
	double uL = UL[1] / rhoL;                   // velocity components
	double vL = UL[2] / rhoL;
	double unL = uL * nx + vL * ny;             // normal velocity
	double qL = sqrt(UL[1] * UL[1] + UL[2] * UL[2]) / rhoL; // velocity magnitude
	double pL = gmi * (UL[3] - 0.5 * rhoL * qL * qL);       // pressure

	// check for non-physical states
	if (pL < 0 || rhoL < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pL << "\n";
		cout << "Density:" << rhoL << "\n";
		exit(1);
	}

	double rHL = UL[3] + pL;   // total enthalpy
	double HL = rHL / rhoL;    // specific enthalpy

	// left flux
	double FL0 = rhoL * unL;
	double FL1 = UL[1] * unL + pL * nx;
	double FL2 = UL[2] * unL + pL * ny;
	double FL3 = rHL * unL;

	// process right state
	double rhoR = UR[0];
	double uR = UR[1] / rhoR;
	double vR = UR[2] / rhoR;
	double unR = uR * nx + vR * ny;
	double qR = sqrt(UR[1] * UR[1] + UR[2] * UR[2]) / rhoR;
	double pR = gmi * (UR[3] - 0.5 * rhoR * qR * qR);

	if (pR < 0 || rhoR < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pR << "\n";
		cout << "Density:" << rhoR << "\n";
		exit(1);
	}

	double rHR = UR[3] + pR;
	double HR = rHR / rhoR;

	// right flux
	double FR0 = rhoR * unR;
	double FR1 = UR[1] * unR + pR * nx;
	double FR2 = UR[2] * unR + pR * ny;
	double FR3 = rHR * unR;

	// difference in states
	double du0 = UR[0] - UL[0], du1 = UR[1] - UL[1];
	double du2 = UR[2] - UL[2], du3 = UR[3] - UL[3];

	// Roe average state
	double di = sqrt(rhoR / rhoL);
	double ui = (di * uR + uL) / (1.0 + di); // roe average velocity components
	double vi = (di * vR + vL) / (1.0 + di);
	double Hi = (di * HR + HL) / (1.0 + di); // roe average specific enthalpy

	double af = 0.5 * (ui * ui + vi * vi);   // roe average kinetic energy (1/2*q^2)
	double ucp = ui * nx + vi * ny;          // roe average normal velocity
	double c2 = gmi * (Hi - af);             // roe average speed of sound squared
	if (c2 < 0) { cout << "Non-physical state!" << "\n"; exit(1); }
	double ci = sqrt(c2);
	double ci1 = 1.0 / ci;

	// eigenvalues
	double l0 = ucp + ci, l1 = ucp - ci, l2 = ucp;

	// entropy fix: parameter = 0.1 * Roe average speed of sound
	double epsilon = ci * .1;
	if (l0 < epsilon && l0 > -epsilon) l0 = 0.5 * (epsilon + l0 * l0 / epsilon);
	if (l1 < epsilon && l1 > -epsilon) l1 = 0.5 * (epsilon + l1 * l1 / epsilon);
	if (l2 < epsilon && l2 > -epsilon) l2 = 0.5 * (epsilon + l2 * l2 / epsilon);

	l0 = fabs(l0); l1 = fabs(l1); l2 = fabs(l2);
	double l3 = l2;

	// average and half-difference of 1st and 2nd eigs
	double s1 = 0.5 * (l0 + l1);
	double s2 = 0.5 * (l0 - l1);

	// left eigenvector product generators
	double G1 = gmi * (af * du0 - ui * du1 - vi * du2 + du3);
	double G2 = -ucp * du0 + du1 * nx + du2 * ny;

	// required functions of G1 and G2 (again, see Theory guide)
	double C1 = G1 * (s1 - l3) * ci1 * ci1 + G2 * s2 * ci1;
	double C2 = G1 * s2 * ci1 + G2 * (s1 - l3);

	// flux assembly
	F[0] = 0.5 * (FL0 + FR0) - 0.5 * (l3 * du0 + C1);
	F[1] = 0.5 * (FL1 + FR1) - 0.5 * (l3 * du1 + C1 * ui + C2 * nx);
	F[2] = 0.5 * (FL2 + FR2) - 0.5 * (l3 * du2 + C1 * vi + C2 * ny);
	F[3] = 0.5 * (FL3 + FR3) - 0.5 * (l3 * du3 + C1 * Hi + C2 * ucp);

	// max wave speed: the eigenvalues are already absolute at this point
	double s_mag = l0;
	if (l1 > s_mag) s_mag = l1;
	if (l2 > s_mag) s_mag = l2;
	return s_mag;
}

inline double rusanovCore(const double* UL, const double* UR, double gamma, double nx, double ny, double* F) {
	// process left state
	double rL = UL[0];                          // density
	double uL = UL[1] / rL;                     // velocity components
	double vL = UL[2] / rL;
	double unL = uL * nx + vL * ny;             // normal velocity
	double qL = sqrt(UL[1] * UL[1] + UL[2] * UL[2]) / rL;   // velocity magnitude
	double pL = (gamma - 1) * (UL[3] - 0.5 * rL * qL * qL); // pressure

	// check for non-physical states
	if (pL < 0 || rL < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pL << "\n";
		cout << "Density:" << rL << "\n";
		exit(1);
	}

	double rHL = UL[3] + pL;              // total enthalpy
	double cL = sqrt(gamma * pL / rL);    // speed of sound

	// left flux
	double FL0 = rL * unL;
	double FL1 = UL[1] * unL + pL * nx;
	double FL2 = UL[2] * unL + pL * ny;
	double FL3 = rHL * unL;

	// process right state
	double rR = UR[0];
	double uR = UR[1] / rR;
	double vR = UR[2] / rR;
	double unR = uR * nx + vR * ny;
	double qR = sqrt(UR[1] * UR[1] + UR[2] * UR[2]) / rR;
	double pR = (gamma - 1) * (UR[3] - 0.5 * rR * qR * qR);

	if (pR < 0 || rR < 0) {
		cout << "Non-physical state!" << "\n";
		exit(1);
	}

	double rHR = UR[3] + pR;
	double cR = sqrt(gamma * pR / rR);

	// right flux
	double FR0 = rR * unR;
	double FR1 = UR[1] * unR + pR * nx;
	double FR2 = UR[2] * unR + pR * ny;
	double FR3 = rHR * unR;

	// max wave speed: |un| + c on each side, so the magnitude does not depend
	// on which way the face normal points (required for conservation)
	double sLmax = fabs(unL) + cL;
	double sRmax = fabs(unR) + cR;
	double s_mag = (sLmax > sRmax) ? sLmax : sRmax;

	// flux assembly
	F[0] = .5 * (FL0 + FR0) - 0.5 * s_mag * (UR[0] - UL[0]);
	F[1] = .5 * (FL1 + FR1) - 0.5 * s_mag * (UR[1] - UL[1]);
	F[2] = .5 * (FL2 + FR2) - 0.5 * s_mag * (UR[2] - UL[2]);
	F[3] = .5 * (FL3 + FR3) - 0.5 * s_mag * (UR[3] - UL[3]);
	return s_mag;
}

inline double hlleCore(const double* UL, const double* UR, double gamma, double nx, double ny, double* F) {
	// process left state
	double rL = UL[0];
	double uL = UL[1] / rL;
	double vL = UL[2] / rL;
	double unL = uL * nx + vL * ny;
	double qL = sqrt(UL[1] * UL[1] + UL[2] * UL[2]) / rL;
	double pL = (gamma - 1) * (UL[3] - 0.5 * rL * qL * qL);

	// check for non-physical states
	if (pL < 0 || rL < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pL << "\n";
		cout << "Density:" << rL << "\n";
		exit(1);
	}

	double rHL = UL[3] + pL;
	double cL = sqrt(gamma * pL / rL);

	// left flux
	double FL0 = rL * unL;
	double FL1 = UL[1] * unL + pL * nx;
	double FL2 = UL[2] * unL + pL * ny;
	double FL3 = rHL * unL;

	// process right state
	double rR = UR[0];
	double uR = UR[1] / rR;
	double vR = UR[2] / rR;
	double unR = uR * nx + vR * ny;
	double qR = sqrt(UR[1] * UR[1] + UR[2] * UR[2]) / rR;
	double pR = (gamma - 1) * (UR[3] - 0.5 * rR * qR * qR);

	if (pR < 0 || rR < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pR << "\n";
		cout << "Density:" << rR << "\n";
		exit(1);
	}

	double rHR = UR[3] + pR;
	double cR = sqrt(gamma * pR / rR);

	// right flux
	double FR0 = rR * unR;
	double FR1 = UR[1] * unR + pR * nx;
	double FR2 = UR[2] * unR + pR * ny;
	double FR3 = rHR * unR;

	// max wave speed
	double smax, smin;
	double sLmax = unL + cL;
	double sRmax = unR + cR;
	double sLmin = unL - cL;
	double sRmin = unR - cR;

	// smax
	if (sLmax < 0 && sRmax < 0) { smax = 0; }
	else { smax = (sLmax > sRmax) ? sLmax : sRmax; }

	// smin
	if (sLmin > 0 && sRmin > 0) { smin = 0; }
	else { smin = (sLmin > sRmin) ? sRmin : sLmin; }

	// largest wave speed magnitude in either direction; smax alone is zero
	// when the flow is entirely along -n, which would allow an infinite dt
	double s_mag = (fabs(smin) > smax) ? fabs(smin) : smax;

	// flux assembly
	double a = 0.5 * ((smax + smin) / (smax - smin));
	double b = (smax * smin) / (smax - smin);
	F[0] = .5 * (FL0 + FR0) - a * (FR0 - FL0) + b * (UR[0] - UL[0]);
	F[1] = .5 * (FL1 + FR1) - a * (FR1 - FL1) + b * (UR[1] - UL[1]);
	F[2] = .5 * (FL2 + FR2) - a * (FR2 - FL2) + b * (UR[2] - UL[2]);
	F[3] = .5 * (FL3 + FR3) - a * (FR3 - FL3) + b * (UR[3] - UL[3]);
	return s_mag;
}

inline double wallFluxCore(const double* u, double nx, double ny, double gam, double* F) {
	double rho = u[0];
	double vx = u[1] / rho, vy = u[2] / rho;
	double vn = vx * nx + vy * ny;
	double vbx = vx - vn * nx;
	double vby = vy - vn * ny;
	double pb = (gam - 1) * (u[3] - .5 * rho * (vbx * vbx + vby * vby));
	if ((pb < 0) || (rho < 0)) { cout << "Non-physical state on the walls!" << "\n"; exit(1); }
	F[0] = 0; F[1] = pb * nx; F[2] = pb * ny; F[3] = 0;
	return sqrt(gam * pb / rho);
}

inline structFlux roe(vector<double> const& UL, vector<double> const& UR, double gamma, vector<double> const& n) {
	// wrapper: allocates a structFlux, so prefer the Core above in hot loops
	structFlux output;
	output.F.resize(4);
	output.s_mag = roeCore(UL.data(), UR.data(), gamma, n[0], n[1], output.F.data());
	return output;
}

/* This function calculates the 2D Rusanov Flux
Inputs: UL, UR = left, right states, as 4x1 vectors
		gamma = ratio of specific heats (e.g. 1.4)
		n = left-to-right unit normal 2x1 vector
Outputs: F = numerical normal flux (4x1 vector)
		 smag = max wave speed estimate
*/
inline structFlux rusanov(vector<double> const& UL, vector<double> const& UR, double gamma, vector<double> const& n) {
	// wrapper: allocates a structFlux, so prefer the Core above in hot loops
	structFlux output;
	output.F.resize(4);
	output.s_mag = rusanovCore(UL.data(), UR.data(), gamma, n[0], n[1], output.F.data());
	return output;
}

/* This function calculates the 2D HLLE Flux
Inputs: UL, UR = left, right states, as 4x1 vectors
		gamma = ratio of specific heats (e.g. 1.4)
		n = left-to-right unit normal 2x1 vector
Outputs: F = numerical normal flux (4x1 vector)
		 smag = max wave speed estimate
*/
inline structFlux HLLE(vector<double> const& UL, vector<double> const& UR, double gamma, vector<double> const& n) {
	// wrapper: allocates a structFlux, so prefer the Core above in hot loops
	structFlux output;
	output.F.resize(4);
	output.s_mag = hlleCore(UL.data(), UR.data(), gamma, n[0], n[1], output.F.data());
	return output;
}

inline structFlux wallFlux(const vector<double>& u, const vector<double>& n, const double& gam) {
	// wrapper: allocates a structFlux, so prefer the Core above in hot loops
	structFlux output;
	output.F.resize(4);
	output.s_mag = wallFluxCore(u.data(), n[0], n[1], gam, output.F.data());
	return output;
}
