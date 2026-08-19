//
//  RK2_FVM.h
//  p2-AEROSP-623
//
//  Created by Jake Yeaman on 2/17/23.
//

#ifndef RK2_FVM_h
#define RK2_FVM_h

#include <iostream>
#include <vector>
#include <cmath>
#include "solver.h"

using namespace std;

/* SSP-RK2 time stepping with local time stepping.
 *
 *     u^(1)   = u^n + dt * f(u^n)
 *     u^(n+1) = u^n + dt/2 * ( f(u^n) + f(u^(1)) )
 *
 * i.e. Heun's method, the optimal two-stage SSP scheme, where f = -R/A. Both
 * stages share the dt built from the stage-1 wave speeds, as a single step must.
 *
 * Non-physical states are not checked here: the flux functions in fluxes.h
 * already report and exit on negative pressure or density.
 */
void rk2(meshData const& mesh, int opt, vector<vector<double>>& u, double Minf,
		 double alphaDeg, string limiterType, double convergedVal, double CFL,
		 int maxIter = 100000) {

	// Hoisted out of the loop: these are overwritten in full every step, so
	// reallocating nelem-by-4 three times per iteration buys nothing.
	vector<vector<double>> f0(mesh.nelem, vector<double>(4, 0.0));
	vector<vector<double>> f1(mesh.nelem, vector<double>(4, 0.0));
	vector<vector<double>> u_f0(mesh.nelem, vector<double>(4, 0.0));
	vector<double> dt(mesh.nelem, 0.0);

	double resL1 = DBL_MAX;
	int niter = 0;
	bool converged = false;

	while (niter < maxIter) {
		niter++;

		// Stage 1. resL1 is the norm of R(u^n), the residual at the current state -
		// this is the quantity proj.pdf Eq. (2) asks to be monitored, and the one
		// FVM1st reports. Testing the intermediate stage's residual instead would
		// converge to a different number.
		vector<vector<double>> residual = secondOrderFV(mesh, opt, u, Minf, alphaDeg, limiterType);
		resL1 = computeL1ResidualNorm(residual);

		if (resL1 <= convergedVal) { converged = true; break; }

		for (int i = 0; i < mesh.nelem; i++) {
			if (residual[i][4] <= 0) {
				cerr << "ERROR: element " << i << " accumulated no wave speed at step "
					 << niter << " - the time step would be infinite\n";
				exit(EXIT_FAILURE);
			}
			dt[i] = (2 * mesh.area[i] * CFL) / residual[i][4]; // local time step

			for (int j = 0; j < 4; j++) {
				f0[i][j]   = -residual[i][j] / mesh.area[i];
				u_f0[i][j] = u[i][j] + dt[i] * f0[i][j];
			}
		}

		// Stage 2, evaluated at the intermediate state and reusing the same dt
		vector<vector<double>> residual2 = secondOrderFV(mesh, opt, u_f0, Minf, alphaDeg, limiterType);

		for (int i = 0; i < mesh.nelem; i++) {
			for (int j = 0; j < 4; j++) {
				f1[i][j] = -residual2[i][j] / mesh.area[i];
				u[i][j]  = u[i][j] + 0.5 * dt[i] * (f0[i][j] + f1[i][j]);
			}
		}

		if (niter % 5 == 0) {
			cout << "Iteration " << niter << " residual is " << resL1 << "\n";
		}
	}

	if (converged) {
		cout << "SSP-RK2 converged in " << niter << " iterations (R = " << resL1 << ")\n";
	} else {
		cout << "Maximum time steps reached (R = " << resL1 << ")\n";
	}

}

#endif /* RK2_FVM_h */


