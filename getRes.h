#include <iostream>
#include <vector>
#include <cmath>
#include "fluxes.h"
using namespace std;

/*function overview*/
vector<vector<double>> getRes(const meshData& mesh, int& opt, vector<vector<double>>& u, 
                              vector<double>& uinf) {
    vector<vector<double>> residual(mesh.nelem, vector<double>(5));
    vector<double> F(4), n(2);
    double s;
 
    // loop over interior edges
    for (int i = 0; i < mesh.interiorFaces.size(); i++){
        int elemL = mesh.interiorFaces[i][2] - 1; // set element left of edge
        int elemR = mesh.interiorFaces[i][3] - 1; // set element right of edge

        // initialize left/right states
        vector<double> uL(4);
        vector<double> uR(4);
        
        // get normal vector for edge
        n[0] = mesh.In[i][0];
        n[1] = mesh.In[i][1];

        // set left/right states
        for (int j = 0; j < 4; j++){
            uL[j] = u[elemL][j]; // set left state
            uR[j] = u[elemR][j]; // set right state
        }
        
        // calculate length of edge
        int node1 = mesh.interiorFaces[i][0] - 1;
        int node2 = mesh.interiorFaces[i][1] - 1;
        double length = sqrt(pow(mesh.nodes[node1][0] - mesh.nodes[node2][0], 2) 
                        + pow(mesh.nodes[node1][1] - mesh.nodes[node2][1], 2));

        // call chosen flux function to compute flux and wave speed using
        // left state, right state, gamma (1.4), and normal vector
        structFlux output;
        if (opt == 1){
            output = roe(uL,uR,1.4,n);
        } else if (opt == 2) {
            output = rusanov(uL,uR,1.4,n);
        } else if (opt == 3) {
            output = HLLE(uL,uR,1.4,n);
        }
        
        F = output.F;
        s = output.s_mag;

        // add F * length to residual of left element, subtract it from right element
        for (int j = 0; j < 4; j++) {
            residual[elemL][j] += F[j] * length;
            residual[elemR][j] -= F[j] * length;
        }
        // // add to the wave speed of left/right elements
        residual[elemL][4] += s * length;
        residual[elemR][4] += s * length;
    }

    // Loop over boundary edges
    for (int i = 0; i < mesh.B2E.size(); i++){
        int elem = mesh.B2E[i][0] - 1; // set element connected to boundary
        vector<double> uTemp(4); // initialize temporary holder for state vector

        // get normal vector for edge
        n[0] = mesh.Bn[i][0];
        n[1] = mesh.Bn[i][1];

        // set temporary holder to element state
        for (int j = 0; j < 4; j++){
            uTemp[j] = u[elem][j];
        } 
        
        // calculate length of edge
        int node1 = mesh.bounds[i][0] - 1;
        int node2 = mesh.bounds[i][1] - 1;
        double length = sqrt(pow(mesh.nodes[node1][0]-mesh.nodes[node2][0],2) 
                             + pow(mesh.nodes[node1][1]-mesh.nodes[node2][1],2));
        
        double gamma = 1.4;
        structFlux output;
        if (mesh.B2E[i][3] == 0) { // 0 = farfield, 1 = wall (set from the .gri group title)

            // Bn points OUT of the domain, i.e. out of uTemp's element, so the
            // interior state is the LEFT one and the freestream is the ghost on
            // the right - the same orientation the interior loop above uses, and
            // the same one the second-order residual uses. Passing (uinf, uTemp) put
            // the exterior state on the left while leaving the normal pointing
            // outward, which negates the central part of the flux while leaving
            // the dissipation alone (F(a,b,n) = 0.5(f(a)+f(b)).n - 0.5s(b-a), so
            // swapping a and b is not a sign flip).
            if (opt == 1) {
                output = roe(uTemp,uinf,gamma,n);
            }
            else if (opt == 2){
                output = rusanov(uTemp,uinf,gamma,n);
            }
            else if (opt == 3){
                output = HLLE(uTemp,uinf,gamma,n);
            }
        } else {
            output = wallFlux(uTemp,n,gamma);
        }
        F = output.F;
        s = output.s_mag;

        // add F*length to residual. R_i = sum over faces of Fhat.n_out * dl, and
        // n_out points out of this element on a boundary face, so this ADDS -
        // exactly as the interior loop adds to its left element. Subtracting here
        // reversed the wall pressure force and the farfield flux.
        for (int j = 0; j < 4; j++) {
            residual[elem][j] += F[j] * length;
        }
        // add wave speed
        residual[elem][4] += s * length;
    }

    return residual;
}
