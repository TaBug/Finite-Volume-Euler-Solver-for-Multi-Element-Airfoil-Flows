#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include "fluxes.h"
#include "solver.h"
using namespace std;

/* Test the Roe, Rusanov, and HLLE fluxes over a boundary with a 30 degree angle of attack */
void unitTest(const vector<double> &uInf, double gamma = 1.4) {
	vector<double> U1 = uInf; // left state
	vector<double> n = {0.5, sqrt(3) / 2};

	cout << "Running unit tests for flux functions with freestream state: [";
	for (double val : uInf) {
		cout << val << " ";
	}
	cout << "]\n";

	// Roe Flux
	structFlux roeFlux = roe(U1, U1, gamma, n);
	cout << "Roe Flux: [";
	for (double val : roeFlux.F) {
		cout << val << " ";
	}
	cout << "]\n";

	// Rusanov flux
	structFlux rusanovflux = rusanov(U1, U1, gamma, n);
	cout << "Rusanov Flux: [";
	for (double val : rusanovflux.F) {
		cout << val << " ";
	}
	cout << "]\n";

	// HLLE Flux
	structFlux hlleFlux = HLLE(U1, U1, gamma, n);
	cout << "HLLE Flux: [";
	for (double val : hlleFlux.F) {
		cout << val << " ";
	}
	cout << "]\n";
}

int main() {
	double gamma = 1.4; // ratio of specific heats
	double Minf = 0.8; // freestream Mach number
	double alpha = 8; // angle of attack in radians

	vector<double> uInf = computeFreestreamState(Minf, alpha); // freestream state vector

	unitTest(uInf, gamma); // run unit tests for flux functions
}
