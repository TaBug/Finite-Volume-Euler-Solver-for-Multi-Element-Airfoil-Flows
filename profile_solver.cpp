// Profiling harness: times the phases of the 2nd-order solver.  Not part of the
// solver build - compile it on its own with .vscode/build.bat.
#include <chrono>
#include <cstdio>
#include "solver.h"
#include "processMesh.h"
#include "data_conversion.h"

using Clock = chrono::high_resolution_clock;
static double ms(Clock::time_point a, Clock::time_point b) {
    return chrono::duration<double, milli>(b - a).count();
}

// Instrumented clone of secondOrderFV in solver.h.  Kept line-for-line identical
// so the split below reflects the real driver.
struct Phase { double geom1=0, accum1=0, areaDiv=0, geom2=0, centroid=0, recon=0, bstate=0, flux=0, scatter=0, alloc=0; };

vector<vector<double>> secondOrderFV_timed(meshData const& mesh, int opt, vector<vector<double>>& u,
                                           double Minf, double alphaDeg, Phase& p) {
    int nelem = mesh.nelem;
    int nfaces = mesh.niedge + mesh.nbedge;
    auto t0 = Clock::now();
    vector<vector<Vector2d>> grad_u(nelem, vector<Vector2d>(4, Vector2d::Zero()));
    p.alloc += ms(t0, Clock::now());

    for (int iFace = 0; iFace < nfaces; iFace++) {
        auto a = Clock::now();
        faceGeom f = computeFaceGeom(mesh, iFace);
        auto b = Clock::now(); p.geom1 += ms(a, b);

        Vector2d nOutL(f.n[0], f.n[1]);
        vector<double> const& uL = u[f.iElemL];
        vector<double> uR = f.isBound
            ? computeBoundaryState(uL, mesh.bounds[f.iFaceLocal][3] == 1, Minf, alphaDeg, mesh.Bn, f.iFaceLocal)
            : u[f.iElemR];
        for (int iU = 0; iU < 4; iU++) {
            Vector2d contribution = 0.5 * (uL[iU] + uR[iU]) * f.length * nOutL;
            grad_u[f.iElemL][iU] += contribution;
            if (!f.isBound) grad_u[f.iElemR][iU] -= contribution;
        }
        p.accum1 += ms(b, Clock::now());
    }

    auto t1 = Clock::now();
    for (int iElem = 0; iElem < nelem; iElem++)
        for (int iU = 0; iU < 4; iU++) grad_u[iElem][iU] /= mesh.area[iElem];
    p.areaDiv += ms(t1, Clock::now());

    auto t2 = Clock::now();
    vector<vector<double>> residuals(nelem, vector<double>(5, 0));
    p.alloc += ms(t2, Clock::now());

    for (int iFace = 0; iFace < nfaces; iFace++) {
        auto a = Clock::now();
        faceGeom f = computeFaceGeom(mesh, iFace);
        auto b = Clock::now(); p.geom2 += ms(a, b);

        Vector2d centroidL = elementCentroid(mesh.nodes, mesh.elem, f.iElemL);
        auto c = Clock::now(); p.centroid += ms(b, c);

        vector<double> uL = u[f.iElemL];
        for (int iU = 0; iU < 4; iU++) uL[iU] += grad_u[f.iElemL][iU].dot(f.midpoint - centroidL);
        vector<double> uR;
        auto d = c;
        if (!f.isBound) {
            Vector2d centroidR = elementCentroid(mesh.nodes, mesh.elem, f.iElemR);
            uR = u[f.iElemR];
            for (int iU = 0; iU < 4; iU++) uR[iU] += grad_u[f.iElemR][iU].dot(f.midpoint - centroidR);
            d = Clock::now(); p.recon += ms(c, d);
        } else {
            d = Clock::now(); p.recon += ms(c, d);
            uR = computeBoundaryState(uL, mesh.bounds[f.iFaceLocal][3] == 1, Minf, alphaDeg, mesh.Bn, f.iFaceLocal);
            auto e2 = Clock::now(); p.bstate += ms(d, e2); d = e2;
        }

        structFlux output = computeFlux(opt, uL, uR, f.n);
        auto e = Clock::now(); p.flux += ms(d, e);

        for (int j = 0; j < 4; j++) residuals[f.iElemL][j] += output.F[j] * f.length;
        residuals[f.iElemL][4] += output.s_mag * f.length;
        if (!f.isBound) {
            for (int j = 0; j < 4; j++) residuals[f.iElemR][j] -= output.F[j] * f.length;
            residuals[f.iElemR][4] += output.s_mag * f.length;
        }
        p.scatter += ms(e, Clock::now());
    }
    return residuals;
}

int main(int argc, char** argv) {
    string griName = (argc > 1) ? argv[1] : "smoothed_local_all.gri";
    int nIter = (argc > 2) ? atoi(argv[2]) : 200;
    int opt = 2;                    // rusanov, matching dat/
    double Minf = 0.25, alphaDeg = 8.0, CFL = 0.9;

    auto t0 = Clock::now();
    string fileNameGri = findGriFile(griName);
    meshData mesh;
    readGriFile(fileNameGri, mesh);
    auto t1 = Clock::now();
    buildMeshTopology(mesh);
    auto t2 = Clock::now();

    printf("mesh: %s  nnode=%d nelem=%d niedge=%d nbedge=%d\n",
           griName.c_str(), mesh.nnode, mesh.nelem, mesh.niedge, mesh.nbedge);
    printf("readGriFile        : %10.2f ms\n", ms(t0, t1));
    printf("buildMeshTopology  : %10.2f ms\n", ms(t1, t2));

    vector<double> uinf = computeFreestreamState(Minf, alphaDeg);
    vector<vector<double>> u(mesh.nelem, uinf);

    // one uninstrumented secondOrderFV, for the timer overhead comparison
    auto t3 = Clock::now();
    for (int k = 0; k < nIter; k++) {
        vector<vector<double>> r = secondOrderFV(mesh, opt, u, Minf, alphaDeg, "NONE");
        (void)r[0][0];
    }
    auto t4 = Clock::now();
    printf("\nsecondOrderFV x%d  : %10.2f ms  (%.4f ms/call, untimed)\n",
           nIter, ms(t3, t4), ms(t3, t4) / nIter);

    // RK2 step cost: two residual evaluations plus the update loops
    auto t5 = Clock::now();
    vector<vector<double>> f0(mesh.nelem, vector<double>(4, 0.0));
    vector<vector<double>> u_f0(mesh.nelem, vector<double>(4, 0.0));
    vector<double> dt(mesh.nelem, 0.0);
    double updateMs = 0, resMs = 0, normMs = 0;
    for (int k = 0; k < nIter; k++) {
        auto a = Clock::now();
        vector<vector<double>> residual = secondOrderFV(mesh, opt, u, Minf, alphaDeg, "NONE");
        auto b = Clock::now(); resMs += ms(a, b);
        double resL1 = computeL1ResidualNorm(residual);
        (void)resL1;
        auto c = Clock::now(); normMs += ms(b, c);
        for (int i = 0; i < mesh.nelem; i++) {
            dt[i] = (2 * mesh.area[i] * CFL) / residual[i][4];
            for (int j = 0; j < 4; j++) {
                f0[i][j] = -residual[i][j] / mesh.area[i];
                u_f0[i][j] = u[i][j] + dt[i] * f0[i][j];
            }
        }
        auto d = Clock::now(); updateMs += ms(c, d);
        vector<vector<double>> residual2 = secondOrderFV(mesh, opt, u_f0, Minf, alphaDeg, "NONE");
        auto e = Clock::now(); resMs += ms(d, e);
        for (int i = 0; i < mesh.nelem; i++)
            for (int j = 0; j < 4; j++)
                u[i][j] = u[i][j] + 0.5 * dt[i] * (-residual2[i][j] / mesh.area[i] + f0[i][j]);
        updateMs += ms(e, Clock::now());
    }
    auto t6 = Clock::now();
    printf("\n--- rk2, %d steps: %.2f ms total (%.4f ms/step) ---\n", nIter, ms(t5, t6), ms(t5, t6)/nIter);
    printf("  secondOrderFV (2 per step) : %10.2f ms  %5.1f%%\n", resMs, 100*resMs/ms(t5,t6));
    printf("  computeL1ResidualNorm      : %10.2f ms  %5.1f%%\n", normMs, 100*normMs/ms(t5,t6));
    printf("  state update loops         : %10.2f ms  %5.1f%%\n", updateMs, 100*updateMs/ms(t5,t6));

    // Instrumented breakdown inside one residual evaluation
    Phase p;
    vector<vector<double>> u2(mesh.nelem, uinf);
    auto t7 = Clock::now();
    for (int k = 0; k < nIter; k++) {
        vector<vector<double>> r = secondOrderFV_timed(mesh, opt, u2, Minf, alphaDeg, p);
        (void)r[0][0];
    }
    auto t8 = Clock::now();
    double tot = ms(t7, t8);
    printf("\n--- inside secondOrderFV, %d calls: %.2f ms wall (timer overhead included) ---\n", nIter, tot);
    printf("  pass1 computeFaceGeom      : %10.2f ms  %5.1f%%\n", p.geom1,   100*p.geom1/tot);
    printf("  pass1 gradient accumulate  : %10.2f ms  %5.1f%%\n", p.accum1,  100*p.accum1/tot);
    printf("  divide gradient by area    : %10.2f ms  %5.1f%%\n", p.areaDiv, 100*p.areaDiv/tot);
    printf("  pass2 computeFaceGeom      : %10.2f ms  %5.1f%%\n", p.geom2,   100*p.geom2/tot);
    printf("  pass2 elementCentroid (L)  : %10.2f ms  %5.1f%%\n", p.centroid,100*p.centroid/tot);
    printf("  pass2 reconstruction       : %10.2f ms  %5.1f%%\n", p.recon,   100*p.recon/tot);
    printf("  pass2 computeBoundaryState : %10.2f ms  %5.1f%%\n", p.bstate,  100*p.bstate/tot);
    printf("  pass2 computeFlux          : %10.2f ms  %5.1f%%\n", p.flux,    100*p.flux/tot);
    printf("  pass2 residual scatter     : %10.2f ms  %5.1f%%\n", p.scatter, 100*p.scatter/tot);
    printf("  container allocation       : %10.2f ms  %5.1f%%\n", p.alloc,   100*p.alloc/tot);
    return 0;
}
