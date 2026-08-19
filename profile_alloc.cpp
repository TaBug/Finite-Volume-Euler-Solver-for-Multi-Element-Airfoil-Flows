// Counts heap allocations per secondOrderFV call, to separate "arithmetic cost"
// from "small-vector churn".  Standalone; not part of the solver build.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <new>

static size_t g_nalloc = 0, g_bytes = 0;
static bool g_count = false;

void* operator new(size_t n) {
    if (g_count) { g_nalloc++; g_bytes += n; }
    void* p = malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    return p;
}
void operator delete(void* p) noexcept { free(p); }
void operator delete(void* p, size_t) noexcept { free(p); }

#include "solver.h"
#include "processMesh.h"

using Clock = chrono::high_resolution_clock;
static double ms(Clock::time_point a, Clock::time_point b) {
    return chrono::duration<double, milli>(b - a).count();
}

int main(int argc, char** argv) {
    string griName = (argc > 1) ? argv[1] : "smoothed_local_all.gri";
    int nIter = (argc > 2) ? atoi(argv[2]) : 100;
    int opt = 2;
    double Minf = 0.25, alphaDeg = 8.0;

    meshData mesh;
    readGriFile(findGriFile(griName), mesh);
    auto tb0 = Clock::now();
    buildMeshTopology(mesh);
    auto tb1 = Clock::now();

    vector<double> uinf = computeFreestreamState(Minf, alphaDeg);
    vector<vector<double>> u(mesh.nelem, uinf);
    int nfaces = mesh.niedge + mesh.nbedge;

    // warm up so the timing excludes first-touch page faults
    { vector<vector<double>> r = secondOrderFV(mesh, opt, u, Minf, alphaDeg, "NONE"); (void)r[0][0]; }

    g_count = true;
    auto t0 = Clock::now();
    for (int k = 0; k < nIter; k++) {
        vector<vector<double>> r = secondOrderFV(mesh, opt, u, Minf, alphaDeg, "NONE");
        (void)r[0][0];
    }
    auto t1 = Clock::now();
    g_count = false;
    size_t refAlloc = g_nalloc, refBytes = g_bytes;

    // the optimized path: geometry and workspace built once, outside the count
    meshGeom geom = buildMeshGeom(mesh, Minf, alphaDeg);
    fvWorkspace ws(mesh.nelem);
    vector<double> uFlat(size_t(mesh.nelem) * 4), R(size_t(mesh.nelem) * 5);
    for (int e = 0; e < mesh.nelem; e++)
        for (int j = 0; j < 4; j++) uFlat[size_t(e) * 4 + j] = u[e][j];
    secondOrderResidual(mesh, geom, ws, opt, uFlat.data(), R.data(), "NONE");  // warm up

    g_nalloc = 0; g_bytes = 0;
    g_count = true;
    auto t2 = Clock::now();
    for (int k = 0; k < nIter; k++)
        secondOrderResidual(mesh, geom, ws, opt, uFlat.data(), R.data(), "NONE");
    auto t3 = Clock::now();
    g_count = false;
    size_t fastAlloc = g_nalloc, fastBytes = g_bytes;

    printf("mesh %-24s nelem=%6d nfaces=%6d   buildMeshTopology %8.1f ms\n",
           griName.c_str(), mesh.nelem, nfaces, ms(tb0, tb1));
    printf("  reference  secondOrderFV       : %8.4f ms/call  %9.0f allocations/call (%.1f per face), %.0f bytes\n",
           ms(t0, t1) / nIter, double(refAlloc) / nIter, double(refAlloc) / nIter / nfaces, double(refBytes) / nIter);
    printf("  optimized  secondOrderResidual : %8.4f ms/call  %9.0f allocations/call, %.0f bytes\n",
           ms(t2, t3) / nIter, double(fastAlloc) / nIter, double(fastBytes) / nIter);
    printf("  speedup: %.1fx\n", ms(t0, t1) / ms(t2, t3));
    return 0;
}
