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

inline structFlux roe(vector<double>& UL, vector<double>& UR, double gamma, vector<double> n) {
	/* INPUTS: UL, UR = left, right states, as 4x1 vectors
			   gamma = ratio of specific heats (e.g. 1.4)
	           n = left-to-right unit normal 2x1 vector
	   OUTPUTS: F = numerical normal flux (4x1 vector)
	            s_mag = max wave speed estimate
	*/
	const double gmi = gamma - 1.0;

	// process left state variables
	double rhoL, uL, vL, unL, qL, pL, rHL, HL;
	rhoL = UL[0]; // density
	uL = UL[1] / rhoL; // velocity components
	vL = UL[2] / rhoL; // velocity components
	unL = uL * n[0] + vL * n[1]; // normal velocity
	qL = sqrt(pow(uL, 2) + pow(vL, 2)); // velocity magnitude
	pL = gmi * (UL[3] - 0.5 * rhoL * pow(qL, 2)); // pressure

	// check for non-physical states
	if (pL < 0 || rhoL < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pL << "\n";
		cout << "Density:" << rhoL << "\n";
	}

	// ensure positive pressure and density
	if (pL < 0) {pL = -pL;}
	if (rhoL < 0) {rhoL = -rhoL;}

	rHL = UL[3] + pL; // total enthalpy
	HL = rHL / rhoL; // specific enthalpy

	// left flux
	vector<double> FL(4, 0.0);
	FL[0] = rhoL * unL;
	FL[1] = UL[1] * unL + pL * n[0];
	FL[2] = UL[2] * unL + pL * n[1];
	FL[3] = rHL * unL;

	// process right state
	double rhoR, uR, vR, unR, qR, pR, rHR, HR;
	rhoR = UR[0]; // density
	uR = UR[1] / rhoR; // velocity components
	vR = UR[2] / rhoR; // velocity components 
	unR = uR * n[0] + vR * n[1]; // normal velocity
	qR = sqrt(pow(uR, 2) + pow(vR, 2)); // velocity magnitude
	pR = gmi * (UR[3] - 0.5 * rhoR * pow(qR, 2)); // pressure

	// check for non-physical states
	if (pR < 0 || rhoR < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pR << "\n";
		cout << "Density:" << rhoR << "\n";
		exit(1);
	}

	// ensure positive pressure and density
	if (pR < 0) {pR = -pR;}
	if (rhoR < 0) {rhoR = -rhoR;}

	rHR = UR[3] + pR; // total enthalpy
	HR = rHR / rhoR; // specific enthalpy

	// right flux
	vector<double> FR(4, 0.0);
	FR[0] = rhoR * unR;
	FR[1] = UR[1] * unR + pR * n[0];
	FR[2] = UR[2] * unR + pR * n[1];
	FR[3] = rHR * unR;

	// difference in states
	vector<double> du;
	du = subtractVectors(UR, UL);

	// Roe average state
	double di, ui, vi, Hi, af, ucp, c2, ci, ci1;
	di = sqrt(rhoR / rhoL);
	ui = (di * uR + uL) / (1.0 + di); // roe average velocity components
	vi = (di * vR + vL) / (1.0 + di); // roe average velocity components
	Hi = (di * HR + HL) / (1.0 + di); // roe average specific enthalpy

	af = 0.5 * (ui * ui + vi * vi); // roe average kinetic energy (1/2*q^2)
	ucp = ui * n[0] + vi * n[1]; // roe average normal velocity
	c2 = gmi * (Hi - af); // roe average speed of sound squared
	if (c2 < 0) {cout << "Non-physical state!" << "\n"; c2 = -c2; exit(1);}
	// if (c2 < 1e-14) {c2 = 1e-14;} // keeps ci1 = 1/ci finite in a near-vacuum state
	ci = sqrt(c2);
	ci1 = 1.0 / ci;

	// eigenvalues
	vector<double> l(3, 0.0);
	l[0] = ucp + ci; l[1] = ucp - ci; l[2] = ucp;

	// entropy fix
	double epsilon, l3;
	epsilon = ci * .1; // entropy fix parameter = 0.1 * Roe average speed of sound
	for (int i = 0; i < 3; i++) {
		if ((l[i] < epsilon) && (l[i] > -epsilon)) {
			l[i] = 0.5 * (epsilon + l[i] * l[i] / epsilon);
		}
	}

	l = absolute(l); l3 = l[2];

	double s1, s2, G1, G2, C1, C2;
	// average and half-difference of 1st and 2nd eigs
	s1 = 0.5 * (l[0] + l[1]);
	s2 = 0.5 * (l[0] - l[1]);

	// left eigenvector product generators
	G1 = gmi * (af * du[0] - ui * du[1] - vi * du[2] + du[3]);
	G2 = -ucp * du[0] + du[1] * n[0] + du[2] * n[1];

	// required functions of G1 and G2(again, see Theory guide)
	C1 = G1 * (s1 - l3) * ci1 * ci1 + G2 * s2 * ci1;
	C2 = G1 * s2 * ci1 + G2 * (s1 - l3);

	// flux assembly
	vector<double> F(4, 0.0);
	F[0] = 0.5 * (FL[0] + FR[0]) - 0.5 * (l3 * du[0] + C1);
	F[1] = 0.5 * (FL[1] + FR[1]) - 0.5 * (l3 * du[1] + C1 * ui + C2 * n[0]);
	F[2] = 0.5 * (FL[2] + FR[2]) - 0.5 * (l3 * du[2] + C1 * vi + C2 * n[1]);
	F[3] = 0.5 * (FL[3] + FR[3]) - 0.5 * (l3 * du[3] + C1 * Hi + C2 * ucp);

	// max wave speed
	double s_mag;
	s_mag = max(l);

	structFlux output;
	output.F = F;
	output.s_mag = s_mag;
	return output;
}

/* This function calculates the 2D Rusanov Flux
Inputs: UL, UR = left, right states, as 4x1 vectors
		gamma = ratio of specific heats (e.g. 1.4)
		n = left-to-right unit normal 2x1 vector
Outputs: F = numerical normal flux (4x1 vector)
		 smag = max wave speed estimate
*/
inline structFlux rusanov(vector<double>& UL, vector<double>& UR, double gamma, vector<double>& n) {
	// process left state
	double rL, uL, vL, unL, qL, pL, rHL, cL;

	rL = UL[0]; // density
	uL = UL[1] / rL; // velocity components
	vL = UL[2] / rL; // velocity components
	unL = uL * n[0] + vL * n[1]; // normal velocity
	qL = sqrt(pow(UL[1], 2) + pow(UL[2], 2)) / rL; // velocity magnitude
	pL = (gamma - 1) * (UL[3] - 0.5 * rL * pow(qL, 2)); // pressure

	// check for non-physical states
	if (pL < 0 || rL < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pL << "\n";
		cout << "Density:" << rL << "\n";
		exit(1);
	}

	// ensure positive pressure and density
	if (pL < 0) {pL = -pL;}
	if (rL < 0) {rL = -rL;}

	rHL = UL[3] + pL; // total enthalpy
	cL = sqrt(gamma * pL / rL); // speed of sound

	// left flux 
	vector<double> FL(4, 0.0);
	FL[0] = rL * unL;
	FL[1] = UL[1] * unL + pL * n[0];
	FL[2] = UL[2] * unL + pL * n[1];
	FL[3] = rHL * unL;

	// process right state
	double rR, uR, vR, unR, qR, pR, rHR, cR;

	rR = UR[0];
	uR = UR[1] / rR;
	vR = UR[2] / rR;
	unR = uR * n[0] + vR * n[1];
	qR = sqrt(pow(UR[1], 2) + pow(UR[2], 2)) / rR;
	pR = (gamma - 1) * (UR[3] - 0.5 * rR * pow(qR, 2));
	if (pR < 0 || rR < 0) {
		cout << "Non-physical state!" << "\n";
		exit(1);
	}
	// if (pR < 0) {
	// 	pR = -pR;
	// }
	// if (rR < 0) {
	// 	rR = -rR;
	// }
	rHR = UR[3] + pR;
	cR = sqrt(gamma * pR / rR);

	// right flux
	vector<double> FR(4, 0.0);
	FR[0] = rR * unR;
	FR[1] = UR[1] * unR + pR * n[0];
	FR[2] = UR[2] * unR + pR * n[1];
	FR[3] = rHR * unR;

	// max wave speed: |un| + c on each side, so the magnitude does not depend
	// on which way the face normal points (required for conservation)
	double s_mag;
	double sLmax = abs(uL * n[0] + vL * n[1]) + cL;
	double sRmax = abs(uR * n[0] + vR * n[1]) + cR;
	s_mag = (sLmax > sRmax) ? sLmax : sRmax;

	// flux assembly
	vector<double> F(4, 0.0);
	F[0] = .5 * (FL[0] + FR[0]) - 0.5 * s_mag * (UR[0] - UL[0]);
	F[1] = .5 * (FL[1] + FR[1]) - 0.5 * s_mag * (UR[1] - UL[1]);
	F[2] = .5 * (FL[2] + FR[2]) - 0.5 * s_mag * (UR[2] - UL[2]);
	F[3] = .5 * (FL[3] + FR[3]) - 0.5 * s_mag * (UR[3] - UL[3]);

	structFlux output;
	output.F = F;
	output.s_mag = s_mag;
	return output;
}

/* This function calculates the 2D HLLE Flux
Inputs: UL, UR = left, right states, as 4x1 vectors
		gamma = ratio of specific heats (e.g. 1.4)
		n = left-to-right unit normal 2x1 vector
Outputs: F = numerical normal flux (4x1 vector)
		 smag = max wave speed estimate
*/
inline structFlux HLLE(vector<double>& UL, vector<double>& UR, double gamma, vector<double>& n) {
	// process left state
	double rL, uL, vL, unL, qL, pL, rHL, cL;

	rL = UL[0];
	uL = UL[1] / rL;
	vL = UL[2] / rL;
	unL = uL * n[0] + vL * n[1];
	qL = sqrt(pow(UL[1], 2) + pow(UL[2], 2)) / rL;
	pL = (gamma - 1) * (UL[3] - 0.5 * rL * pow(qL, 2));

	// check for non-physical states
	if (pL < 0 || rL < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pL << "\n";
		cout << "Density:" << rL << "\n";
		exit(1);
	}
	// if (pL < 0) {
	// 	pL = -pL;
	// }
	// if (rL < 0) {
	// 	rL = -rL;
	// }
	rHL = UL[3] + pL;
	cL = sqrt(gamma * pL / rL);

	// left flux 
	vector<double> FL(4, 0.0);
	FL[0] = rL * unL;
	FL[1] = UL[1] * unL + pL * n[0];
	FL[2] = UL[2] * unL + pL * n[1];
	FL[3] = rHL * unL;

	// process right state
	double rR, uR, vR, unR, qR, pR, rHR, cR;

	rR = UR[0];
	uR = UR[1] / rR;
	vR = UR[2] / rR;
	unR = uR * n[0] + vR * n[1];
	qR = sqrt(pow(UR[1], 2) + pow(UR[2], 2)) / rR;
	pR = (gamma - 1) * (UR[3] - 0.5 * rR * pow(qR, 2));
	if (pR < 0 || rR < 0) {
		cout << "Non-physical state!" << "\n";
		cout << "Pressure:" << pR << "\n";
		cout << "Density:" << rR << "\n";
		exit(1);
	}
	// if (pR < 0) {
	// 	pR = -pR;
	// }
	// if (rR < 0) {
	// 	rR = -rR;
	// }
	rHR = UR[3] + pR;
	cR = sqrt(gamma * pR / rR);

	// right flux
	vector<double> FR(4, 0.0);
	FR[0] = rR * unR;
	FR[1] = UR[1] * unR + pR * n[0];
	FR[2] = UR[2] * unR + pR * n[1];
	FR[3] = rHR * unR;

	// max wave speed
	double s_mag, smax, smin;
	double sLmax = uL * n[0] + vL * n[1] + cL;
	double sRmax = uR * n[0] + vR * n[1] + cR;
	double sLmin = uL * n[0] + vL * n[1] - cL;
	double sRmin = uR * n[0] + vR * n[1] - cR;

	//smax
	if (sLmax < 0 && sRmax < 0) {smax=0;} 
	else { 
		if (sLmax > sRmax) {smax = sLmax;}
	else {smax = sRmax;}}

	// smin
	if (sLmin > 0 && sRmin > 0){smin=0;}
	else {
		if (sLmin > sRmin) {smin = sRmin;}
	else {smin = sLmin;}}	
	
	// largest wave speed magnitude in either direction; smax alone is zero
	// when the flow is entirely along -n, which would allow an infinite dt
	s_mag = (abs(smin) > smax) ? abs(smin) : smax;
	
	// flux assembly
	vector<double> F(4, 0.0);
	F[0] = .5 * (FL[0] + FR[0]) - 0.5 * ((smax+smin)/(smax-smin)) * (FR[0] - FL[0])+ ((smax*smin)/(smax-smin)) * (UR[0] - UL[0]);
	F[1] = .5 * (FL[1] + FR[1]) - 0.5 * ((smax+smin)/(smax-smin)) * (FR[1] - FL[1])+ ((smax*smin)/(smax-smin)) * (UR[1] - UL[1]);
	F[2] = .5 * (FL[2] + FR[2]) - 0.5 * ((smax+smin)/(smax-smin)) * (FR[2] - FL[2])+ ((smax*smin)/(smax-smin)) * (UR[2] - UL[2]);
	F[3] = .5 * (FL[3] + FR[3]) - 0.5 * ((smax+smin)/(smax-smin)) * (FR[3] - FL[3])+ ((smax*smin)/(smax-smin)) * (UR[3] - UL[3]);

	structFlux output;
	output.F = F;
	output.s_mag = s_mag;
	return output;
}

inline structFlux wallFlux(const vector<double>& u, const vector<double>& n, const double& gam) {
	vector<double> v(2);
    vector<double> vb(2);
   	vector<double> F(4);
    double pb, s_mag;

	// local copy of the density: u is const, and the sign fix below has to be
	// applied before it reaches the sqrt in the wave speed
	double rho = u[0];

	v = { u[1] / rho, u[2] / rho };

	vb[0] = v[0] - (v[0] * n[0] + v[1] * n[1]) * n[0];
	vb[1] = v[1] - (v[0] * n[0] + v[1] * n[1]) * n[1];
	pb = (gam - 1) * (u[3] - .5 * rho * (pow(vb[0], 2) + pow(vb[1], 2)));
	if ((pb < 0) || (rho < 0)) {cout << "Non-physical state on the walls!" << "\n"; exit(1);}
	// if (pb < 0) { pb = -pb; }
	// if (rho < 0) { rho = -rho; }
	F = { 0,pb * n[0],pb * n[1],0 };
	s_mag = sqrt(gam * pb / rho);
    
	// return the flux and the wave speed
    structFlux output;
	output.F = F;
	output.s_mag = s_mag;
	return output;
}
