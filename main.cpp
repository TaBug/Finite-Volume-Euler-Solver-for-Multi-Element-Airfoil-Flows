#include <filesystem>
#include "solver.h"
#include "processMesh.h"
#include "FE_FVM.h"
#include "RK2_FVM.h"
#include "elem2Edge.h"
using namespace std;

int main() {
    // provided condition
    double gamma = 1.4; // ratio of specific heats
    double Minf = 0.25; // freestream Mach number
    double alphaDeg = 8.0; // angle of attack in degrees
    vector<double> uinf = computeFreestreamState(Minf, alphaDeg); // freestream state vector [rho, rho*u, rho*v, rho*E]
 
    // find the .gri file in the working directory
    string fileNameGri = findGriFile("smoothed_local_all.gri"); 
    cout << "Mesh File (.gri): " << fileNameGri << "\n";
    
    // read the .gri file and construct the mesh data structure
    meshData mesh;
    readGriFile(fileNameGri, mesh); // nodes, elements and boundaries
    buildMeshTopology(mesh); // interiorFaces, I2E, B2E, In, Bn, Area

    // initialize the state to free-stream condition
    vector<vector<double>> u (mesh.nelem, computeFreestreamState(Minf, alphaDeg));

    int opt;     // flux function choice, read from the user by FVM1st
    double CFL;  // likewise, and reused by rk2 below
    FVM1st(mesh.bounds, mesh.nodes, mesh.interiorFaces, u, mesh.B2E, mesh.Bn, mesh.In, mesh.nelem, opt, uinf, CFL);
    
    vector<vector<double>> U_firstOrder = u;
    
    vector<vector<int>> elemBounds = genElemBounds(me1sh.elem, mesh.I2E, mesh.B2E);
    vector<vector<int>> globalEdges = genGlobalEdge(mesh.bounds, mesh.interiorFaces);
    
    cout << opt;
    rk2(opt, U_firstOrder, mesh.Area, mesh.nodes, mesh.elem, Minf, alphaDeg, mesh.Bn, mesh.In, elemBounds, mesh.bounds, mesh.interiorFaces, globalEdges, mesh.I2E, mesh.B2E, "BJ", pow(10,-2), mesh.Area, mesh.elem.size(), CFL);

    return 0;
}
