#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include "fluxes.h"
#include "getRes.h"
using namespace std;

// First-order finite volume solve, marched to steady state with Forward Euler
// and local time stepping. Returns when the L1 residual drops below resTol, or
// when the user's iteration cap is reached.
inline void FVM1st(meshData& mesh, vector<vector<double>>& u, int& opt, vector<double> uinf, double& CFL, const double resTol) {
    // CFL is an out-parameter so the value entered here also reaches rk2 in
    // main; it used to be a local, leaving main's copy uninitialized
    int maxIter;

    // User input for which Flux Function
    cout << "Choose Flux Function:\n" << "1 - Roe\n" << "2 - Rusanov\n" << "3 - HLLE\n" << "Enter Option: ";
    cin >> opt;

    // User input for CFL value
    cout << "Enter CFL: ";
    cin >> CFL;

    // User input for maximum time iterations
    cout << "Enter Max. Time Iterations: ";
    cin >> maxIter;

    // A failed extraction leaves its target at 0, which would otherwise surface
    // much later as a solve that never moves (CFL = 0) or one that stops on the
    // very first step (maxIter = 0)
    if (!cin || opt < 1 || opt > 3 || CFL <= 0 || maxIter <= 0) {
        cerr << "ERROR: invalid input (opt=" << opt << ", CFL=" << CFL
             << ", maxIter=" << maxIter << ")\n";
        exit(EXIT_FAILURE);
    }

    cout << "Running first-order finite volume method...\n";
    // Time stepping until the residual is smaller than the given tolerance
    double resL1 = DBL_MAX;
    int t = 0;
    while (resL1 > resTol) {
        resL1 = 0;

        // calculate residual of each element
        vector<vector<double>> residual = getRes(mesh, opt, u, uinf);

        // update state using Forward Euler (first order accurate). Note that
        // dtOverArea is dt/A, not dt: the cell area in dt = 2*CFL*A/sum(s*L)
        // cancels against the 1/A in dU/dt = -R/A, so it never appears here.
        for (int i = 0; i < mesh.nelem; i++) {
            if (residual[i][4] <= 0) {
                cerr << "ERROR: element " << i << " accumulated no wave speed at step "
                     << t << " - the time step would be infinite\n";
                exit(EXIT_FAILURE);
            }
            double dtOverArea = (2 * CFL) / residual[i][4]; // local time step / area
            for (int j = 0; j < 4; j++) {
                u[i][j] = u[i][j] - (dtOverArea * residual[i][j]); // update state
            }
        }

        // calculate L1 Residual
        for (int i = 0; i < mesh.nelem; i++) {
            for (int j = 0; j < 4; j++) {
                resL1 += abs(residual[i][j]);
            }
        }

        // NaN compares false against everything, so without this check the loop
        // condition above would fall straight through and report a diverged run
        // as having converged
        if (!isfinite(resL1)) {
            cerr << "ERROR: residual is " << resL1 << " at step " << t
                 << " - the solve has diverged\n";
            exit(EXIT_FAILURE);
        }

        t++;
        cout << resL1 << "\n";

        if (t >= maxIter) {
            cout << "Maximum time steps reached (R = " << resL1 << ")\n";
            return; // hand control back to the caller instead of killing the process
        }
    }

    cout << "First-order finite volume method has converged! (time step = " << t << ") \n";
}
