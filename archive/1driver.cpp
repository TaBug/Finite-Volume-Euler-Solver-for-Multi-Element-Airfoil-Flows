#include <iostream>
#include <vector>
#include <cmath>
#include "FE_FVM.h"
#include <iostream>
#include "processMesh (2).h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <numeric>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

int main () {
    cout << "Computational Fluid Dynamics II Driver...\n";
    int nnode, nbedge, nelem;
    string fileName = "/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/CoarseMesh (1)";
    getMatSizes(fileName, nnode,nbedge,nelem);
    
    int niedge = 0.5*(3*nelem - nbedge);
    
    vector<vector<double>> nodes(nnode, vector<double>(2));
    vector<vector<double>> elem(nelem, vector<double>(3)); // [node1, node2, node3]
    vector<vector<double>> bounds(nbedge, vector<double>(3)); // [node1, node2, boundaryGroup]
    
    readGmshFile(fileName,nnode,nbedge,nelem,nodes,elem,bounds);
    generateGri("/mnt/c/Users/sasha/Documents/MATLAB/aero623/.vscode/Project2-main/CoarseMesh.gri", nnode, nbedge, nelem, nodes, elem, bounds);
    vector<vector<double>> interiorFaces = genInteriorFaceVec(niedge, nelem, nbedge, bounds, elem);
    
    niedge = int(interiorFaces.size());

    // GENERATING TASK 2 DATA STRUCTURES
    vector<vector<double>> I2E = genI2E(niedge, nelem, interiorFaces, elem);
    vector<vector<double>> B2E = genB2E(nbedge, nelem, bounds, elem);
    vector<vector<double>> In = genIn(niedge, interiorFaces, elem, nodes, I2E);
    vector<vector<double>> Bn = genBn(nbedge, bounds, nodes, elem, B2E);
    vector<double> Area = genArea(nelem, elem, nodes);

    vector<vector<double>> u(nelem, vector<double> (4));
    for(int i=0;i<nelem;++i)  {
        u[i][0]=1;
        u[i][1]=.25*cos(8*3.14159/180);
        u[i][2]=.25*sin(8*3.14159/180);
        u[i][3]=1/(.4*1.4)+.25*.25/2;
    }

    FVM_1st(bounds,nodes,interiorFaces,u,B2E,Bn,In,nelem,Area);
    ofstream outputFile;
    outputFile.open("output.txt");
    
    for(int iElem = 0; iElem < u.size(); iElem++){
        for(int iU = 0; iU < 1; iU++){
            outputFile << u[iElem][iU] << " ";
        }
        outputFile << "\n";
    }
    outputFile.close();
    return 0;
}