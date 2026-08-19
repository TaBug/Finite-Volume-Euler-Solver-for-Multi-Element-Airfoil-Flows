// Control experiment: the same 2nd-order residual with the per-face heap traffic
// removed (face geometry cached once, states in std::array, flux written in
// place).  Verifies bit-for-bit-close agreement with secondOrderFV, then times
// both, so the cost attributed to allocation in the phase breakdown can be
// checked against an actual speedup.  Standalone; not part of the solver build.
#include <chrono>
#include <cstdio>
#include <array>
#include <cmath>
#include "solver.h"
#include "processMesh.h"

using Clock = chrono::high_resolution_clock;
static double ms(Clock::time_point a, Clock::time_point b) {
    return chrono::duration<double, milli>(b - a).count();
}

struct FaceCache {
    int iElemL, iElemR;   // iElemR = -1 on a boundary face
    bool isWall;
    int iFaceLocal;
    double nx, ny, len;
    double mx, my;        // face midpoint
};

static vector<FaceCache> buildFaceCache(meshData const& mesh) {
    int nfaces = mesh.niedge + mesh.nbedge;
    vector<FaceCache> fc(nfaces);
    for (int i = 0; i < nfaces; i++) {
        faceGeom f = computeFaceGeom(mesh, i);
        fc[i] = { f.iElemL, f.iElemR,
                  f.isBound && mesh.bounds[f.iFaceLocal][3] == 1,
                  f.iFaceLocal, f.n[0], f.n[1], f.length,
                  f.midpoint[0], f.midpoint[1] };
    }
    return fc;
}

// Centroids do not move, so they are cached with the face geometry.
static vector<array<double,2>> buildCentroids(meshData const& mesh) {
    vector<array<double,2>> c(mesh.nelem);
    for (int e = 0; e < mesh.nelem; e++) {
        Vector2d v = elementCentroid(mesh.nodes, mesh.elem, e);
        c[e] = { v[0], v[1] };
    }
    return c;
}

static inline double rusanovFast(const double* UL, const double* UR, double nx, double ny, double* F) {
    const double g = 1.4;
    double rL = UL[0], uL = UL[1]/rL, vL = UL[2]/rL;
    double unL = uL*nx + vL*ny;
    double qL = sqrt(UL[1]*UL[1] + UL[2]*UL[2]) / rL;
    double pL = (g-1)*(UL[3] - 0.5*rL*qL*qL);
    double rHL = UL[3] + pL, cL = sqrt(g*pL/rL);
    double FL0 = rL*unL, FL1 = UL[1]*unL + pL*nx, FL2 = UL[2]*unL + pL*ny, FL3 = rHL*unL;

    double rR = UR[0], uR = UR[1]/rR, vR = UR[2]/rR;
    double unR = uR*nx + vR*ny;
    double qR = sqrt(UR[1]*UR[1] + UR[2]*UR[2]) / rR;
    double pR = (g-1)*(UR[3] - 0.5*rR*qR*qR);
    double rHR = UR[3] + pR, cR = sqrt(g*pR/rR);
    double FR0 = rR*unR, FR1 = UR[1]*unR + pR*nx, FR2 = UR[2]*unR + pR*ny, FR3 = rHR*unR;

    double sL = fabs(unL) + cL, sR = fabs(unR) + cR;
    double s = (sL > sR) ? sL : sR;
    F[0] = .5*(FL0+FR0) - 0.5*s*(UR[0]-UL[0]);
    F[1] = .5*(FL1+FR1) - 0.5*s*(UR[1]-UL[1]);
    F[2] = .5*(FL2+FR2) - 0.5*s*(UR[2]-UL[2]);
    F[3] = .5*(FL3+FR3) - 0.5*s*(UR[3]-UL[3]);
    return s;
}

static inline void boundaryStateFast(const double* Uc, bool isWall, const double* uinf,
                                     double nx, double ny, double* Ub) {
    if (!isWall) { for (int i = 0; i < 4; i++) Ub[i] = uinf[i]; return; }
    double rho = Uc[0], vx = Uc[1]/rho, vy = Uc[2]/rho;
    double vn = vx*nx + vy*ny;
    Ub[0] = rho; Ub[1] = rho*(vx - vn*nx); Ub[2] = rho*(vy - vn*ny); Ub[3] = Uc[3];
}

// flat [nelem*4] state, flat [nelem*5] residual - no per-face heap traffic
static void secondOrderFV_fast(meshData const& mesh, vector<FaceCache> const& fc,
                               vector<array<double,2>> const& cent, const double* uinf,
                               const double* u, double* R, vector<double>& grad) {
    int nelem = mesh.nelem;
    fill(grad.begin(), grad.end(), 0.0);          // [nelem][4][2]
    fill(R, R + size_t(nelem)*5, 0.0);

    for (size_t i = 0; i < fc.size(); i++) {
        FaceCache const& f = fc[i];
        const double* uL = u + size_t(f.iElemL)*4;
        double ubuf[4];
        const double* uR;
        if (f.iElemR < 0) { boundaryStateFast(uL, f.isWall, uinf, f.nx, f.ny, ubuf); uR = ubuf; }
        else               uR = u + size_t(f.iElemR)*4;

        for (int k = 0; k < 4; k++) {
            double c = 0.5 * (uL[k] + uR[k]) * f.len;
            grad[size_t(f.iElemL)*8 + k*2 + 0] += c * f.nx;
            grad[size_t(f.iElemL)*8 + k*2 + 1] += c * f.ny;
            if (f.iElemR >= 0) {
                grad[size_t(f.iElemR)*8 + k*2 + 0] -= c * f.nx;
                grad[size_t(f.iElemR)*8 + k*2 + 1] -= c * f.ny;
            }
        }
    }
    for (int e = 0; e < nelem; e++) {
        double inv = 1.0 / mesh.area[e];
        for (int k = 0; k < 8; k++) grad[size_t(e)*8 + k] *= inv;
    }

    for (size_t i = 0; i < fc.size(); i++) {
        FaceCache const& f = fc[i];
        double uLr[4], uRr[4];
        double dxL = f.mx - cent[f.iElemL][0], dyL = f.my - cent[f.iElemL][1];
        for (int k = 0; k < 4; k++)
            uLr[k] = u[size_t(f.iElemL)*4 + k]
                   + grad[size_t(f.iElemL)*8 + k*2 + 0]*dxL + grad[size_t(f.iElemL)*8 + k*2 + 1]*dyL;

        if (f.iElemR >= 0) {
            double dxR = f.mx - cent[f.iElemR][0], dyR = f.my - cent[f.iElemR][1];
            for (int k = 0; k < 4; k++)
                uRr[k] = u[size_t(f.iElemR)*4 + k]
                       + grad[size_t(f.iElemR)*8 + k*2 + 0]*dxR + grad[size_t(f.iElemR)*8 + k*2 + 1]*dyR;
        } else {
            boundaryStateFast(uLr, f.isWall, uinf, f.nx, f.ny, uRr);
        }

        double F[4];
        double s = rusanovFast(uLr, uRr, f.nx, f.ny, F);
        for (int k = 0; k < 4; k++) R[size_t(f.iElemL)*5 + k] += F[k] * f.len;
        R[size_t(f.iElemL)*5 + 4] += s * f.len;
        if (f.iElemR >= 0) {
            for (int k = 0; k < 4; k++) R[size_t(f.iElemR)*5 + k] -= F[k] * f.len;
            R[size_t(f.iElemR)*5 + 4] += s * f.len;
        }
    }
}

int main(int argc, char** argv) {
    string griName = (argc > 1) ? argv[1] : "smoothed_local_all.gri";
    int nIter = (argc > 2) ? atoi(argv[2]) : 200;
    int opt = 2;
    double Minf = 0.25, alphaDeg = 8.0;

    meshData mesh;
    readGriFile(findGriFile(griName), mesh);
    buildMeshTopology(mesh);

    vector<double> uinf = computeFreestreamState(Minf, alphaDeg);
    vector<vector<double>> u(mesh.nelem, uinf);
    vector<double> uflat(size_t(mesh.nelem)*4);
    for (int e = 0; e < mesh.nelem; e++)
        for (int k = 0; k < 4; k++) uflat[size_t(e)*4+k] = u[e][k];

    auto tc0 = Clock::now();
    vector<FaceCache> fc = buildFaceCache(mesh);
    vector<array<double,2>> cent = buildCentroids(mesh);
    auto tc1 = Clock::now();

    vector<double> R(size_t(mesh.nelem)*5), grad(size_t(mesh.nelem)*8);

    // correctness check against the shipped driver
    vector<vector<double>> Rref = secondOrderFV(mesh, opt, u, Minf, alphaDeg, "NONE");
    secondOrderFV_fast(mesh, fc, cent, uinf.data(), uflat.data(), R.data(), grad);
    // scaled by the largest residual in the field: individual entries pass through
    // zero, so a per-entry relative error is meaningless there
    double maxabs = 0, scale = 0;
    for (int e = 0; e < mesh.nelem; e++)
        for (int k = 0; k < 5; k++) {
            maxabs = max(maxabs, fabs(Rref[e][k] - R[size_t(e)*5+k]));
            scale  = max(scale,  fabs(Rref[e][k]));
        }
    double maxerr = maxabs / scale;

    auto t0 = Clock::now();
    for (int k = 0; k < nIter; k++) {
        vector<vector<double>> r = secondOrderFV(mesh, opt, u, Minf, alphaDeg, "NONE");
        (void)r[0][0];
    }
    auto t1 = Clock::now();
    for (int k = 0; k < nIter; k++)
        secondOrderFV_fast(mesh, fc, cent, uinf.data(), uflat.data(), R.data(), grad);
    auto t2 = Clock::now();

    double a = ms(t0, t1) / nIter, b = ms(t1, t2) / nIter;
    printf("mesh %-24s nelem=%6d nfaces=%6d\n", griName.c_str(), mesh.nelem, mesh.niedge + mesh.nbedge);
    printf("  max relative residual difference : %.3e\n", maxerr);
    printf("  one-time geometry precompute     : %8.3f ms\n", ms(tc0, tc1));
    printf("  secondOrderFV (as written)       : %8.4f ms/call\n", a);
    printf("  same math, no per-face heap churn: %8.4f ms/call   (%.1fx faster)\n", b, a / b);
    return 0;
}
