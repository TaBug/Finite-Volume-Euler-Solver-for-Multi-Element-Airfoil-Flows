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
 *
 * u is a flat [nelem][4] array (stride 4) and geom is the precomputed face
 * geometry from buildMeshGeom. Every buffer the loop needs is allocated once
 * below: at 100k iterations, allocating them per step was the dominant cost.
 */
void rk2(meshData const& mesh, meshGeom const& geom, int opt, vector<double>& u,
		 string limiterType, double convergedVal, double CFL, int maxIter = 100000) {

	const int nelem = mesh.nelem;

	// Hoisted out of the loop: these are overwritten in full every step, so
	// reallocating them per iteration buys nothing. The two residual buffers and
	// the gradient workspace used to be allocated and freed inside the residual
	// call itself, twice per step.
	fvWorkspace ws(nelem);
	vector<double> residual(size_t(nelem) * 5, 0.0);
	vector<double> residual2(size_t(nelem) * 5, 0.0);
	vector<double> f0(size_t(nelem) * 4, 0.0);
	vector<double> u_f0(size_t(nelem) * 4, 0.0);
	vector<double> dt(nelem, 0.0);

	double resL1 = DBL_MAX;
	int niter = 0;
	bool converged = false;

	while (niter < maxIter) {
		niter++;

		// Stage 1. resL1 is the norm of R(u^n), the residual at the current state -
		// this is the quantity proj.pdf Eq. (2) asks to be monitored, and the one
		// FVM1st reports. Testing the intermediate stage's residual instead would
		// converge to a different number.
		secondOrderResidual(mesh, geom, ws, opt, u.data(), residual.data(), limiterType);
		resL1 = computeL1ResidualNorm(residual.data(), nelem);

		if (resL1 <= convergedVal) { converged = true; break; } // if the residual is small enough, exit the loop

		for (int i = 0; i < nelem; i++) {
			const double* Ri = residual.data() + size_t(i) * 5;
			if (Ri[4] <= 0) {
				cerr << "ERROR: element " << i << " accumulated no wave speed at step "
					 << niter << " - the time step would be infinite\n";
				exit(EXIT_FAILURE);
			}
			dt[i] = (2 * mesh.area[i] * CFL) / Ri[4]; // local time step

			double* f0i = f0.data() + size_t(i) * 4;
			double* ui = u.data() + size_t(i) * 4;
			double* u_f0i = u_f0.data() + size_t(i) * 4;
			for (int j = 0; j < 4; j++) {
				f0i[j]   = -Ri[j] / mesh.area[i];
				u_f0i[j] = ui[j] + dt[i] * f0i[j];
			}
		}

		// Stage 2, evaluated at the intermediate state and reusing the same dt
		secondOrderResidual(mesh, geom, ws, opt, u_f0.data(), residual2.data(), limiterType);

		for (int i = 0; i < nelem; i++) {
			const double* R2i = residual2.data() + size_t(i) * 5;
			const double* f0i = f0.data() + size_t(i) * 4;
			double* ui = u.data() + size_t(i) * 4;
			for (int j = 0; j < 4; j++) {
				double f1 = -R2i[j] / mesh.area[i];
				ui[j] = ui[j] + 0.5 * dt[i] * (f0i[j] + f1);
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
