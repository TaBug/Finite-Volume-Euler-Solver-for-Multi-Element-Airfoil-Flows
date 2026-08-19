#include <filesystem>
#include <sstream>
#include <iomanip>
#include "solver.h"
#include "processMesh.h"
#include "FE_FVM.h"
#include "RK2_FVM.h"
#include "data_conversion.h"
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

    // ask if the user has already run the first-order solver and has a converged solution to use as an initial condition for RK2
    cout << "Have you already run the first-order solver? (y/n): ";
    char firstOrderRun;
    cin >> firstOrderRun;

    int opt;     // flux function choice, read from the user by FVM1st
    double CFL;  // likewise, and reused by rk2 below
    const char* fluxName[] = { "", "roe", "rusanov", "hlle" }; // opt is 1-3, validated below

    // name the files after the flux and CFL so the results are self-describing,
    // and so a reloaded solution can be traced back to the settings that made it.
    // to_string(CFL) is avoided: it always emits 6 decimals, giving "CFL0.700000"
    auto solutionFile = [&](const string &order) {
        ostringstream cflStr;
        cflStr << setprecision(3) << CFL;
        return "dat/" + string(fluxName[opt]) + "_CFL" + cflStr.str() + "_" + order + ".dat";
    };

    // initialize the state to free-stream condition
    vector<vector<double>> u (mesh.nelem, uinf); // [elemID][rho, rho*u, rho*v, rho*E]
    if (firstOrderRun != 'y' && firstOrderRun != 'Y') { // not yet run, so run the first-order solver
        // first-order finite volume method to get a converged solution
        FVM1st(mesh, u, opt, uinf, CFL, 1e-5);

        string filename = solutionFile("firstOrder");
        state2text(u, filename); // save the converged solution to a file
        cout << "First-order solution saved to '" << filename << "'\n";
    } else { // already run
        // read the converged solution from the file to use as an initial condition for RK2
        cout << "Enter the filename of the converged first-order solution: ";
        string filename;
        cin >> filename;
        u = text2state(filename);

        // opt and CFL are set by FVM1st in the branch above, so this branch has to
        // supply them itself - rk2 uses both, and reading them uninitialized is UB.
        // The state file holds only the four conserved variables, so they come from
        // the filename, falling back to the user when it does not match the scheme.
        if (!filename2settings(filename, opt, CFL)) {
            cout << "Could not read the flux and CFL from '" << filename << "'.\n";
            cout << "Choose Flux Function:\n" << "1 - Roe\n" << "2 - Rusanov\n" << "3 - HLLE\n" << "Enter Option: ";
            cin >> opt;
            cout << "Enter CFL: ";
            cin >> CFL;
            if (!cin || opt < 1 || opt > 3 || CFL <= 0) {
                cerr << "ERROR: invalid input (opt=" << opt << ", CFL=" << CFL << ")\n";
                exit(EXIT_FAILURE);
            }
        }
        cout << "Converged first-order solution loaded from '" << filename
             << "' (flux = " << fluxName[opt] << ", CFL = " << CFL << ")\n";
    }

    vector<vector<double>> u_firstOrder = u; // store the converged solution from first-order FVM to use as initial condition for RK2

    cout << "Would you like to proceed with the second-order solver? (y/n): ";
    char proceed;
    cin >> proceed;
    if (proceed != 'y' && proceed != 'Y') {
        cout << "Exiting program.\n";
        return 0;
    }

    cout << "Running second-order Runge-Kutta finite volume method...\n";
    // "NONE" for now: the BJ and LCD limiters are implemented in solver.h but are
    // not yet called from secondOrderFV, which rejects any other value rather than
    // quietly running unlimited.
    rk2(mesh, opt, u_firstOrder, Minf, alphaDeg, "NONE", 1e-5, CFL);
    cout << "Second-order Runge-Kutta finite volume method has converged!\n";

    // rk2 updates u_firstOrder in place and writes nothing itself, so the save
    // happens here - the message below used to announce a file that never existed
    string outFilename = solutionFile("secondOrder");
    state2text(u_firstOrder, outFilename);
    cout << "Solution saved to '" << outFilename << "'\n";
    return 0;
}
