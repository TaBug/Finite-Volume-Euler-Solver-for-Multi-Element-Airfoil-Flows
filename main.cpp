#include <filesystem>
#include <sstream>
#include <iomanip>
#include "solver.h"
#include "processMesh.h"
#include "FE_FVM.h"
#include "RK2_FVM.h"
#include "elem2Edge.h"
#include "data2text.h"
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
    buildMeshTopology(mesh); // interiorFaces, I2E, B2E, In, Bn, area

    // initialize the state to free-stream condition
    vector<vector<double>> u (mesh.nelem, uinf); // [elemID][rho, rho*u, rho*v, rho*E]

    int opt;     // flux function choice, read from the user by FVM1st
    double CFL;  // likewise, and reused by rk2 below
    
    // first-order finite volume method to get a converged solution
    FVM1st(mesh, u, opt, uinf, CFL, 1e-5);

    // name the file after the flux and CFL so the results are self-describing.
    // to_string(CFL) is avoided: it always emits 6 decimals, giving "CFL0.700000"
    const char* fluxName[] = { "", "roe", "rusanov", "hlle" }; // opt is 1-3, validated in FVM1st
    ostringstream cflStr;
    cflStr << setprecision(3) << CFL;
    string filename = "dat/" + string(fluxName[opt]) + "_CFL" + cflStr.str() + "_firstOrder.dat";

    state2text(u, filename); // save the converged solution to a file
    cout << "First-order solution saved to '" << filename << "'\n";

    vector<vector<double>> U_firstOrder = u; // store the converged solution from first-order FVM to use as initial condition for RK2
    vector<vector<int>> elemBounds = genElemBounds(mesh.elem, mesh.I2E, mesh.B2E);
    vector<vector<int>> globalEdges = genGlobalEdge(mesh.bounds, mesh.interiorFaces);
    

    cout << "Would you like to proceed with the second-order solver? (y/n): ";
    char proceed;
    cin >> proceed;
    if (proceed != 'y' && proceed != 'Y') {
        cout << "Exiting program.\n";
        return 0;
    }

    // cout << "Running second-order Runge-Kutta finite volume method...\n";
    // rk2(mesh, opt, U_firstOrder, Minf, alphaDeg, elemBounds, globalEdges, "BJ", pow(10,-2), CFL);
    // cout << "Second-order Runge-Kutta finite volume method has converged!\n";
    // cout << "Solution saved to 'solution.dat'\n";
    // return 0;
}
