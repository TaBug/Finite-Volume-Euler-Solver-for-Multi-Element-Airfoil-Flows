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
    double gamma = 1.4; // ratio of specific heats

    // The two flow conditions the project statement specifies. Both share an
    // angle of attack; only the freestream Mach number differs. Neither is
    // recoverable from a state file afterwards - it holds nothing but the four
    // conserved variables - so the choice is stamped into the solution name and
    // read back from there by scripts/post_process.py.
    cout << "Choose Flow Condition:\n"
         << "1 - Subsonic  (Minf = 0.25, alpha = 8 deg)\n"
         << "2 - Transonic (Minf = 0.50, alpha = 8 deg)\n"
         << "Enter Option: ";
    int flowOption;
    cin >> flowOption;
    if (!cin || (flowOption != 1 && flowOption != 2)) {
        cerr << "ERROR: invalid flow condition (expected 1 = subsonic, 2 = transonic)\n";
        exit(EXIT_FAILURE);
    }
    const double Minf = (flowOption == 1) ? 0.25 : 0.50; // freestream Mach number
    const double alphaDeg = 8.0;                         // angle of attack in degrees
    const string flowName = (flowOption == 1) ? "subsonic" : "transonic";
    cout << "Flow Condition: " << flowName << " (Minf = " << Minf
         << ", alpha = " << alphaDeg << " deg)\n";

    vector<double> uinf = computeFreestreamState(Minf, alphaDeg); // freestream state vector [rho, rho*u, rho*v, rho*E]

    // let the user pick one of the meshes in the gri folder
    string fileNameGri = chooseGriFile();
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

    // The mesh is part of what identifies a run: the same flux and CFL on a
    // different mesh is a different solution, and a state file is only meaningful
    // on the mesh it was written for - hence the element-count check on reload
    // below. stem() drops the "gri/" and the ".gri".
    const string meshName = filesystem::path(fileNameGri).stem().string();

    // name the files after the flux, CFL, flow condition, mesh and - for the
    // second-order runs - the limiter, so the results are self-describing and a
    // reloaded solution can be traced back to the settings that made it:
    //
    //     <flux>_CFL<cfl>_<order>_<flow>_<mesh>[_<limiter>].dat
    //
    // Field order is load-bearing at both ends. filename2settings reads the flux
    // off the front and the CFL after "CFL", so everything past those is free;
    // post_process.py reads the limiter off the LAST underscore-separated field,
    // so mesh and limiter stay at the tail in that order. The flow tag goes
    // between the order and the mesh, which leaves the whole
    // "<flux>_CFL<cfl>_<order>" prefix untouched and keeps the mesh adjacent to
    // the limiter. Mesh names may themselves contain underscores
    // ("smoothed_local_all") - fine, since no parse counts fields from the middle.
    // to_string(CFL) is avoided: it always emits 6 decimals, giving "CFL0.700000"
    auto solutionFile = [&](const string &order, const string &tag = "") {
        ostringstream cflStr;
        cflStr << setprecision(3) << CFL;
        // first- and second-order solutions go to their own subfolders; openForWrite
        // creates the folder, so neither has to exist beforehand
        return "dat/" + order + "/" + string(fluxName[opt]) + "_CFL" + cflStr.str()
               + "_" + order + "_" + flowName + "_" + meshName
               + (tag.empty() ? "" : "_" + tag) + ".dat";
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
        // pick the converged solution to use as an initial condition for RK2
        string filename = chooseDatFile("firstOrder");
        cout << "Solution File (.dat): " << filename << "\n";
        u = text2state(filename);

        // A solution written on a different mesh loads without complaint but leaves
        // every later loop indexing past the end of u. Easy to hit now that the mesh
        // and the solution are picked from two independent menus, so check the shape.
        if (int(u.size()) != mesh.nelem) {
            cerr << "ERROR: '" << filename << "' holds " << u.size()
                 << " elements but the mesh has " << mesh.nelem << "\n"
                 << "       The solution was written on a different mesh.\n";
            exit(EXIT_FAILURE);
        }

        // The flow condition is not recoverable from the state either, so a
        // subsonic restart marched under transonic freestream conditions would
        // run to convergence and quietly produce a wrong answer. The name records
        // which one it was, so say so when it disagrees with what was picked
        // above. Only a warning: a deliberate cross-condition restart is a
        // legitimate way to reach a harder case.
        const string otherFlow = (flowOption == 1) ? "transonic" : "subsonic";
        if (filename.find("_" + otherFlow + "_") != string::npos) {
            cerr << "WARNING: '" << filename << "' was written at the " << otherFlow
                 << " condition,\n         but this run is set up for " << flowName
                 << " (Minf = " << Minf << "). Continuing anyway.\n";
        }

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

    // Plot the first-order solution

    cout << "Would you like to proceed with the second-order solver? (y/n): ";
    char proceed;
    cin >> proceed;
    if (proceed != 'y' && proceed != 'Y') {
        cout << "Exiting program.\n";
        return 0;
    }

    cout << "Running second-order Runge-Kutta finite volume method...\n";

    // Everything about the mesh that the residual needs - face normals, edge
    // lengths, reconstruction offsets, the freestream state - is fixed for the
    // whole solve, so it is built once here rather than per face per iteration.
    meshGeom geom = buildMeshGeom(mesh, Minf, alphaDeg);

    // The second-order path carries the state flat (stride 4); the first-order
    // solver above still uses the nested layout, so convert at the seam.
    vector<double> uFlat = flattenState(u_firstOrder);

    // Choose the limiter type.
    string limiterType = "NONE";
    cout << "Choose Limiter Type:\n" << "1 - None\n" << "2 - Barth-Jespersen\n" << "3 - LCD\n" << "Enter Option: ";
    int limiterOption;
    cin >> limiterOption;
    switch (limiterOption) {
        case 1:
            limiterType = "NONE";
            break;
        case 2:
            limiterType = "BJ";
            break;
        case 3:
            limiterType = "LCD";
            break;
        default:
            cerr << "ERROR: invalid limiter option " << limiterOption << "\n";
            exit(EXIT_FAILURE);
    }

    // User input for maximum time iterations
    int maxIter;
    cout << "Enter Max. Time Iterations: ";
    cin >> maxIter;

    // Run the second-order solver with the converged first-order solution as the initial condition.
    rk2(mesh, geom, opt, uFlat, limiterType, 1e-5, CFL, maxIter);
    cout << "Second-order Runge-Kutta finite volume method has converged!\n";

    // rk2 updates uFlat in place and writes nothing itself, so the save
    // happens here - the message below used to announce a file that never existed
    string outFilename = solutionFile("secondOrder", limiterType);
    state2text(uFlat, outFilename);
    cout << "Solution saved to '" << outFilename << "'\n";

    
    return 0;
}
