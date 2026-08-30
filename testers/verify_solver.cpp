// End-to-end equivalence + timing check for the optimized second-order path.
//
// Marches the SAME initial condition for the same number of SSP-RK2 steps twice:
// once through the reference implementation (secondOrderFV, nested vectors, the
// pre-optimization loop body reproduced verbatim below) and once through the
// shipped fast path (secondOrderResidual + rk2). Reports the largest difference
// in the final state and the time each took.
//
// Standalone; not part of the solver build. Build with build_opt.bat.
#include <chrono>
#include <cstdio>
#include <cmath>
#include "solver.h"
#include "processMesh.h"
#include "RK2_FVM.h"
#include "data_conversion.h"

using Clock = chrono::high_resolution_clock;
static double ms(Clock::time_point a, Clock::time_point b) {
    return chrono::duration<double, milli>(b - a).count();
}

// The original rk2 loop body, kept here as the reference to compare against.
static void rk2_reference(meshData const& mesh, int opt, vector<vector<double>>& u, double Minf,
                          double alphaDeg, string limiterType, double CFL, int nSteps) {
    vector<vector<double>> f0(mesh.nelem, vector<double>(4, 0.0));
    vector<vector<double>> f1(mesh.nelem, vector<double>(4, 0.0));
    vector<vector<double>> u_f0(mesh.nelem, vector<double>(4, 0.0));
    vector<double> dt(mesh.nelem, 0.0);

    for (int niter = 0; niter < nSteps; niter++) {
        vector<vector<double>> residual = secondOrderFV(mesh, opt, u, Minf, alphaDeg, limiterType);
        for (int i = 0; i < mesh.nelem; i++) {
            dt[i] = (2 * mesh.area[i] * CFL) / residual[i][4];
            for (int j = 0; j < 4; j++) {
                f0[i][j]   = -residual[i][j] / mesh.area[i];
                u_f0[i][j] = u[i][j] + dt[i] * f0[i][j];
            }
        }
        vector<vector<double>> residual2 = secondOrderFV(mesh, opt, u_f0, Minf, alphaDeg, limiterType);
        for (int i = 0; i < mesh.nelem; i++) {
            for (int j = 0; j < 4; j++) {
                f1[i][j] = -residual2[i][j] / mesh.area[i];
                u[i][j]  = u[i][j] + 0.5 * dt[i] * (f0[i][j] + f1[i][j]);
            }
        }
    }
}

int main(int argc, char** argv) {
    string griName  = (argc > 1) ? argv[1] : "smoothed_local_all.gri";
    int nSteps      = (argc > 2) ? atoi(argv[2]) : 50;
    string initFile = (argc > 3) ? argv[3] : "";
    double Minf = 0.25, alphaDeg = 8.0, CFL = 0.9;

    meshData mesh;
    readGriFile(findGriFile(griName), mesh);
    buildMeshTopology(mesh);

    // Realistic states matter: starting from uniform freestream would leave every
    // gradient at zero and exercise none of the reconstruction.
    vector<double> uinf = computeFreestreamState(Minf, alphaDeg);
    vector<vector<double>> u0(mesh.nelem, uinf);
    if (!initFile.empty()) {
        vector<vector<double>> loaded = text2state(initFile);
        if (int(loaded.size()) != mesh.nelem) {
            printf("ERROR: %s has %zu rows but the mesh has %d elements\n",
                   initFile.c_str(), loaded.size(), mesh.nelem);
            return 1;
        }
        u0 = loaded;
        printf("initial condition: %s\n", initFile.c_str());
    } else {
        printf("initial condition: uniform freestream\n");
    }

    meshGeom geom = buildMeshGeom(mesh, Minf, alphaDeg);
    const char* fluxName[] = { "", "roe", "rusanov", "hlle" };

    printf("mesh %s  nelem=%d nfaces=%d   %d SSP-RK2 steps, CFL=%.2f\n\n",
           griName.c_str(), mesh.nelem, mesh.niedge + mesh.nbedge, nSteps, CFL);
    printf("%-9s %14s %14s %8s %14s\n", "flux", "reference", "optimized", "speedup", "max |du|/|u|");

    bool allGood = true;
    for (int opt = 1; opt <= 3; opt++) {
        vector<vector<double>> uRef = u0;
        auto t0 = Clock::now();
        rk2_reference(mesh, opt, uRef, Minf, alphaDeg, "NONE", CFL, nSteps);
        auto t1 = Clock::now();

        // convergedVal = 0 so it never exits early; maxIter caps it at nSteps.
        // rk2 prints progress - harmless here, filtered by the caller.
        vector<double> uFast = flattenState(u0);
        auto t2 = Clock::now();
        rk2(mesh, geom, opt, uFast, "NONE", 0.0, CFL, nSteps);
        auto t3 = Clock::now();

        double num = 0, den = 0;
        for (int e = 0; e < mesh.nelem; e++)
            for (int k = 0; k < 4; k++) {
                num = max(num, fabs(uRef[e][k] - uFast[size_t(e) * 4 + k]));
                den = max(den, fabs(uRef[e][k]));
            }
        double rel = num / den;
        if (rel > 1e-10) allGood = false;

        printf("%-9s %11.1f ms %11.1f ms %7.1fx %14.3e\n",
               fluxName[opt], ms(t0, t1), ms(t2, t3), ms(t0, t1) / ms(t2, t3), rel);
    }

    printf("\n%s\n", allGood ? "PASS: optimized path matches the reference to <1e-10"
                             : "FAIL: states diverged");
    return allGood ? 0 : 1;
}
