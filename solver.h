#ifndef solver_h
#define solver_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <numeric>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "processMesh.h"
#include "fluxes.h"

using namespace std;
using namespace Eigen;

// Not M_PI: that is a POSIX extension MSVC gates behind _USE_MATH_DEFINES, which
// only takes effect if defined before the first <cmath> in the translation unit.
// main.cpp includes <filesystem> ahead of this header, so a #define here is dead.
// constexpr double PI = 3.14159265358979323846;

double computeL1ResidualNorm(vector<vector<double>> const& residuals) {
	double L1ResidualNorm = 0; // initialize L1 residual to 0 every time iteration

	// calculate L1 Residual and stop calculation if < 10e-5
	for (int iElem = 0; iElem < residuals.size(); iElem++) {
		for (int iU = 0; iU < 4; iU++) {
			L1ResidualNorm += abs(residuals[iElem][iU]);
		} 
	}

	return L1ResidualNorm;
}

vector<double> computeFreestreamState(const double& Minf, const double& alphaDeg) {
	double alphaRad = alphaDeg * M_PI / 180.0;
	double gamma = 1.4;
	vector<double> uInf = { 1.0, Minf * cos(alphaRad), Minf * sin(alphaRad), (1 / (gamma * gamma - gamma)) + 0.5 * (Minf * Minf) };
	return uInf;
}

vector<Vector3d> computeP(vector<vector<double>> const& nodes, vector<vector<double>> const& elem, vector<double> const& U_cell, int const iCell) {
	double xCentroid = ((nodes[elem[iCell][0] - 1][0]) + (nodes[elem[iCell][1] - 1][0]) + (nodes[elem[iCell][2] - 1][0])) / 3;
	double yCentroid = ((nodes[elem[iCell][0] - 1][1]) + (nodes[elem[iCell][1] - 1][1]) + (nodes[elem[iCell][2] - 1][1])) / 3;
	vector<Vector3d> P;
	P.reserve(U_cell.size());

	for (int i = 0; i < U_cell.size(); i++) {
		P.push_back(Vector3d(xCentroid, yCentroid, U_cell[i]));
	}

	return P;
}

// State on the far side of boundary face iBound, so the interior flux routine can
// be reused unchanged along the domain edge.
//
// Wall: an inviscid wall allows slip but no through-flow, so the normal component
// of the velocity is projected out, v_b = v - (v.n)n. Density and total energy are
// carried over from the interior cell on purpose - holding rho*E fixed while the
// velocity shrinks is what lets the flux routine recover the wall pressure as
// p_b = (gamma-1)(rho*E - 0.5*rho*|v_b|^2).
//
// Bn must hold UNIT normals: the projection only removes the full normal component
// when |n| = 1, and silently under-corrects otherwise.
//
// Farfield: the full freestream state, valid while the outer boundary is far enough
// out that no characteristic information needs to travel back into the domain.
inline vector<double> computeBoundaryState(vector<double> const& U_cell, bool isWall, double Minf, double alphaDeg, vector<vector<double>> const& Bn, int iBound) {

	if (!isWall) { // farfield
		return computeFreestreamState(Minf, alphaDeg);
	}

	double rho = U_cell[0];
	Vector2d vCell = { U_cell[1] / rho, U_cell[2] / rho };
	Vector2d n = { Bn[iBound][0], Bn[iBound][1] };
	Vector2d vb = vCell - vCell.dot(n) * n;

	return { rho, rho * vb[0], rho * vb[1], U_cell[3] };
}

vector<Vector3d> computeP_boundary(vector<vector<double>> const& nodes, vector<vector<double>> const& elem, vector<double> const& U_cell, int const iCell, const int iFace, bool isWall, double Minf, double alphaDeg, vector<vector<double>> const& Bn, int iBound) {

	// Finding the boundary edge midpoint
	vector<size_t> nodeIndices = { 0,1,2 }; // representing localface 1, 2, and 3.
	nodeIndices.erase(nodeIndices.begin() + (iFace - 1)); // removing local face index so our two remaining indices correspond to the index of edge nodes
	vector<double> edgeMidpoint = { (nodes[elem[iCell][nodeIndices[0]] - 1][0] + nodes[elem[iCell][nodeIndices[1]] - 1][0]) / 2, (nodes[elem[iCell][nodeIndices[0]] - 1][1] + nodes[elem[iCell][nodeIndices[1]] - 1][1]) / 2 }; // cell0 and cellk edge interface midpoint

	// Finding current, non-ghost cell center
	double xCentroid = ((nodes[elem[iCell][0] - 1][0]) + (nodes[elem[iCell][1] - 1][0]) + (nodes[elem[iCell][2] - 1][0])) / 3;
	double yCentroid = ((nodes[elem[iCell][0] - 1][1]) + (nodes[elem[iCell][1] - 1][1]) + (nodes[elem[iCell][2] - 1][1])) / 3;

	// Finding ghost cell center
	double dx = edgeMidpoint[0] - xCentroid;
	double dy = edgeMidpoint[1] - yCentroid;
	double xCentroidGhost = edgeMidpoint[0] + dx;
	double yCentroidGhost = edgeMidpoint[1] + dy;

	//    vector<double> U_ghost;
	//    U_ghost.reserve(4);
	//
	//    if(isWall == true){ // Wall Boundary
	//        double rho = U_cell[0];
	//        Vector2d vCell = {U_cell[1]/rho,U_cell[2]/rho};
	//        Vector2d n = {Bn[iBound][0],Bn[iBound][1]};
	//        Vector2d vb = vCell - (vCell.dot(n))*n;
	//        U_ghost = {rho, rho*vb[0], rho*vb[1], U_cell[3]};
	//    }
	//    else{ // Freestream
	//        U_ghost = computeFreestreamState(Minf, alphaDeg);
	//    }

	vector<double> U_ghost = computeBoundaryState(U_cell, isWall, Minf, alphaDeg, Bn, iBound);

	// Computing P
	vector<Vector3d> P;
	P.reserve(4);

	for (int iU = 0; iU < 4; iU++) {
		P.emplace_back(Vector3d(xCentroidGhost, yCentroidGhost, U_ghost[iU]));
	}

	return P;

}

// iNeighbor: a vector of neightboring cell indeces (size of two on boundaries and size of three on interior)
vector<Vector2d> computeL(vector<vector<double>> const& nodes, vector<vector<double>> const& elem, vector<vector<double>> const& U, int iCell, vector<int> const iNeighbor, vector<int> const iFaces, double Minf, double alphaDeg, vector<vector<double>> const& Bn, vector<vector<double>> const& bounds, vector<vector<double>> const& interiorFaces, vector<vector<int>> const& elemBounds) {
	// iFaces: vector consisting of the local face index of the interface edge between cell0 and cellk, same ordering as iNeighbor

	// TODO: acount for boundary edge states: farfield = free-stream and wall = density and energy equivalent to the state inside the domain
	// Assuming that for boundary cases, L_123 are corrected to the only possible L vector value. Therefore, size should always be zero

	vector<Vector3d> Pi;
	vector<Vector3d> Pj;
	vector<Vector3d> Pk;
	vector<Vector3d> zerosVec(4, Vector3d::Zero());


	// LET iNeighbor BE NEGATIVE IFF ON A WALL BOUNDARY
	// TODO: IF iNeighbor is NEGATIVE, corresponding iFace should be the local face of iCell. Ensure this is passed in correctly to the code
	// If not on a wall boundary, compute P with the neighboring cell coordinates and neighboring cell states
	if (iNeighbor[0] < 0) {

		int iGlobal = elemBounds[iCell][iFaces[0] - 1];
		iGlobal2Local currFace = iG2L(iGlobal, bounds, interiorFaces);
		bool isWall = false;
		if (bounds[currFace.index][3] == 1) { // 1 = wall, 0 = farfield (from the .gri group title)
			isWall = true;
		}
		Pi = computeP_boundary(nodes, elem, U[iCell], iCell, iFaces[0], isWall, Minf, alphaDeg, Bn, currFace.index);

	}

	else {
		Pi = computeP(nodes, elem, U[iNeighbor[0]], iNeighbor[0]);
	}

	if (iNeighbor[1] < 0) {
		int iGlobal = elemBounds[iCell][iFaces[1] - 1];
		iGlobal2Local currFace = iG2L(iGlobal, bounds, interiorFaces);
		bool isWall = false;
		if (bounds[currFace.index][3] == 1) { // 1 = wall, 0 = farfield (from the .gri group title)
			isWall = true;
		}
		Pj = computeP_boundary(nodes, elem, U[iCell], iCell, iFaces[1], isWall, Minf, alphaDeg, Bn, currFace.index);
	}

	else {
		Pj = computeP(nodes, elem, U[iNeighbor[1]], iNeighbor[1]);
	}

	if (iNeighbor[2] < 0) {
		int iGlobal = elemBounds[iCell][iFaces[2] - 1];
		iGlobal2Local currFace = iG2L(iGlobal, bounds, interiorFaces);
		bool isWall = false;
		if (bounds[currFace.index][3] == 1) { // 1 = wall, 0 = farfield (from the .gri group title)
			isWall = true;
		}
		Pk = computeP_boundary(nodes, elem, U[iCell], iCell, iFaces[2], isWall, Minf, alphaDeg, Bn, currFace.index);

	}

	else {
		Pk = computeP(nodes, elem, U[iNeighbor[2]], iNeighbor[2]);
	}


	vector<Vector3d> n;
	n.reserve(Pi.size());

	for (int i = 0; i < Pi.size(); i++) {
		n.emplace_back(Vector3d((Pj[i] - Pi[i]).cross(Pk[i] - Pi[i])));
	}

	vector<Vector2d> Lijk;
	Lijk.reserve(Pi.size());

	for (int i = 0; i < Pi.size(); i++) {
		Lijk.emplace_back(Vector2d(-n[i][0] / n[i][2], -n[i][1] / n[i][2]));
	}


	return Lijk;

}

vector<Vector2d> computeL_LCD(vector<Vector2d> L, vector<vector<double>> const& U, vector<vector<double>> const& nodes, vector<vector<double>> const& elem, int iCell, vector<int> const iNeighbor, vector<int> const iFaces, double Minf, double alphaDeg, vector<vector<double>> const& Bn, vector<vector<int>> const& elemBounds, vector<vector<double>> const& bounds, vector<vector<double>> const& interiorFaces) {
	// TODO: IF iNeighbor is NEGATIVE, corresponding iFace should be the local face of iCell. Ensure this is passed in correctly to the code
	// iNeighbor: vector consisting of the index of the neighboring elements
	// iFaces: vector consisting of the local face index of the interface edge between cell0 and cellk, same ordering as iNeighbor

	vector<Vector2d> Llcd;
	Llcd.reserve(U[0].size());

	vector<double> alpha(U[0].size(), DBL_MAX); // storing an independent alpha value for each state variable
	vector<double> cellCentroid = { (nodes[elem[iCell][0] - 1][0] + nodes[elem[iCell][1] - 1][0] + nodes[elem[iCell][2] - 1][0]) / 3,(nodes[elem[iCell][0] - 1][1] + nodes[elem[iCell][1] - 1][1] + nodes[elem[iCell][2] - 1][1]) / 3 };

	for (int i = 0; i < 3; i++) {

		vector<size_t> nodeIndices = { 0,1,2 }; // representing localface 1, 2, and 3.
		nodeIndices.erase(nodeIndices.begin() + (iFaces[i] - 1)); // removing local face index so our two remaining indices correspond to the index of edge nodes

		int ik = iNeighbor[i]; // cellk index
		vector<double> edgeMidpoint;

		if (ik < 0) { // if on a boudnary
			edgeMidpoint = { (nodes[elem[iCell][nodeIndices[0]] - 1][0] + nodes[elem[iCell][nodeIndices[1]] - 1][0]) / 2, (nodes[elem[iCell][nodeIndices[0]] - 1][1] + nodes[elem[iCell][nodeIndices[1]] - 1][1]) / 2 }; // cell0 and boundary edge interface midpoint
		}
		else { // If not on a boundary
			edgeMidpoint = { (nodes[elem[iCell][nodeIndices[0]] - 1][0] + nodes[elem[iCell][nodeIndices[1]] - 1][0]) / 2, (nodes[elem[iCell][nodeIndices[0]] - 1][1] + nodes[elem[iCell][nodeIndices[1]] - 1][1]) / 2 }; // cell0 and cellk edge interface midpoint
		}

		Vector2d r0k(edgeMidpoint[0] - cellCentroid[0], edgeMidpoint[1] - cellCentroid[1]);

		for (int j = 0; j < alpha.size(); j++) { // populating alpha for each state

			double alpha_k;
			double stateDifference;
			double rDotL;

			if (ik < 0) { // if on a boundary
				int iGlobal = elemBounds[iCell][iFaces[i] - 1];
				iGlobal2Local currFace = iG2L(iGlobal, bounds, interiorFaces);
				bool isWall = false;
				if (bounds[currFace.index][3] == 1) { // 1 = wall, 0 = farfield (from the .gri group title)
					isWall = true;
				}
				vector<double> Ub = computeBoundaryState(U[iCell], isWall, Minf, alphaDeg, Bn, currFace.index);
				stateDifference = Ub[j] - U[iCell][j];
				rDotL = r0k.dot(L[j]);
			}
			else { // if not on a boundary
				stateDifference = U[ik][j] - U[iCell][j];
				rDotL = r0k.dot(L[j]);
			}

			if (rDotL > max(stateDifference, 0.0)) {
				alpha_k = max(stateDifference, 0.0) / rDotL;
			}
			else if (r0k.dot(L[j]) < min(stateDifference, 0.0)) {
				alpha_k = min(stateDifference, 0.0) / rDotL;
			}
			else {
				alpha_k = 1;
			}

			// Updating alpha if alphak < alphaCurrent

			if (alpha_k < alpha[j]) {
				alpha[j] = alpha_k;
			}
		}


	}

	Llcd = L;
	for (int i = 0; i < alpha.size(); i++) {
		Llcd[i] *= alpha[i];
	}

	return Llcd;

}

/// BARTH JESPERSON

vector<Vector2d> compute_rN(vector<vector<double>> const& nodes, vector<vector<double>> const& elem, int iCell) {

	vector<double> cellCentroid = { (nodes[elem[iCell][0] - 1][0] + nodes[elem[iCell][1] - 1][0] + nodes[elem[iCell][2] - 1][0]) / 3,(nodes[elem[iCell][0] - 1][1] + nodes[elem[iCell][1] - 1][1] + nodes[elem[iCell][2] - 1][1]) / 3 };

	// Zero(), not Vector2d(2, 0.0): the two-argument form is Eigen's (x, y)
	// constructor, so that would seed every entry with the vector (2, 0).
	vector<Vector2d> rN(3, Vector2d::Zero());

	// Iterate through each node
	for (int i = 0; i < 3; i++) {
		double xN = nodes[elem[iCell][i] - 1][0];
		double yN = nodes[elem[iCell][i] - 1][1];
		rN[i][0] = xN - cellCentroid[0]; // xComponent of rN
		rN[i][1] = yN - cellCentroid[1]; // yComponent of rN
	}

	return rN;

}

vector<Vector2d> barthJespersen(vector<vector<double>> const& nodes, vector<vector<double>> const& elem, vector<double> const& area, vector<vector<double>> const& U, int iCell, vector<int> const iNeighbor, double Minf, double alphaDeg, vector<vector<double>> const& Bn, vector<vector<int>> const& elemBounds, vector<vector<double>> const& bounds, vector<vector<double>> const& interiorFaces, vector<int> const iFaces) {
	// vector<vector<double>> uALL(4, vector<double>(4,0.0)); // EACH ROW REPRESENTS A SPECIFIC STATE VARIABLE (E.G DENSTIY)
	vector<vector<double>> uALL(4, vector<double>(1 + iNeighbor.size(), 0.0));

	// Populating uALL for iCell
	uALL[0][0] = U[iCell][0];
	uALL[1][0] = U[iCell][1];
	uALL[2][0] = U[iCell][2];
	uALL[3][0] = U[iCell][3];

	// Populating uALL for neighbors
	for (int i = 0; i < iNeighbor.size(); i++) {

		if (iNeighbor[i] < 0) { // Boundary
			if (iNeighbor[i] == -1) { // freestream
				vector<double> uFS = computeFreestreamState(Minf, alphaDeg);
				uALL[0][i + 1] = uFS[0];
				uALL[1][i + 1] = uFS[1];
				uALL[2][i + 1] = uFS[2];
				uALL[3][i + 1] = uFS[3];
			}
			else { // wall
				int iFaceGlobal = elemBounds[iCell][iFaces[i] - 1];
				iGlobal2Local iFaceInfo = iG2L(iFaceGlobal, bounds, interiorFaces);
				int iFaceLocal = iFaceInfo.index;
				vector<double> uWall = computeBoundaryState(U[iCell], true, Minf, alphaDeg, Bn, iFaceLocal);
				uALL[0][i + 1] = uWall[0];
				uALL[1][i + 1] = uWall[1];
				uALL[2][i + 1] = uWall[2];
				uALL[3][i + 1] = uWall[3];
			}
		}
		else {
			uALL[0][i + 1] = U[iNeighbor[i]][0];
			uALL[1][i + 1] = U[iNeighbor[i]][1];
			uALL[2][i + 1] = U[iNeighbor[i]][2];
			uALL[3][i + 1] = U[iNeighbor[i]][3];
		}
	}

	vector<double> uMAX(4, 0.0);
	vector<double> uMIN(4, 0.0);

	// Iterate through each state variable
	for (int i = 0; i < 4; i++) {
		uMAX[i] = *(max_element(uALL[i].begin(), uALL[i].end()));
		uMIN[i] = *(min_element(uALL[i].begin(), uALL[i].end()));
	}


	// Find cell0 state
	vector<double> u0 = U[iCell];


	// Obtain all rNs
	vector<Vector2d> rN = compute_rN(nodes, elem, iCell);

	//out rN
	// cout << "rN" << endl;
	// for (int i = 0; i < rN.size(); i++) {
	// // Loop through each element in the inner vector
	// for (int j = 0; j < rN[i].size(); j++) {
	//     cout << rN[i][j] << " ";
	// }
	// cout << endl;
	// }

// Optain all L
	vector<Vector2d> L = computeL(nodes, elem, U, iCell, iNeighbor, iFaces, Minf, alphaDeg, Bn, bounds, interiorFaces, elemBounds);

	//out L
	// cout << "L" << endl;
	// for (int i = 0; i < L.size(); i++) {
	// // Loop through each element in the inner vector
	// for (int j = 0; j < L[i].size(); j++) {
	//     cout << L[i][j] << " ";
	// }
	// cout << endl;
	// }

// Obtain node states
// Each row is a node and each col is a state variable
	vector<vector<double>> uiN;
	uiN.reserve(3);

	// Iterate through each node
	for (int iN = 0; iN < 3; iN++) {
		vector<double> currState;
		currState.reserve(4);
		// Iterate through each state var.
		for (int iU = 0; iU < 4; iU++) {
			currState.emplace_back(u0[iU] + rN[iN].dot(L[iU]));
		}
		uiN.emplace_back(currState);

	}

	//out Uin
	// cout << "Uin" << endl;
	// for (int i = 0; i < uiN.size(); i++) {
	// // Loop through each element in the inner vector
	// for (int j = 0; j < uiN[i].size(); j++) {
	//     cout << uiN[i][j] << " ";
	// }
	// cout << endl;
	// }


// Computing alphaNs
	vector<double> alphas(4, DBL_MAX); // Final, minimul alpha values

	for (int iN = 0; iN < 3; iN++) { // Iterate through each node

		for (int iU = 0; iU < 4; iU++) { // Iterate through each state variable

			if ((uiN[iN][iU] - u0[iU]) > 0) {

				double alpha_iN = min(1.0, (uMAX[iU] - u0[iU]) / (uiN[iN][iU] - u0[iU]));

				if (alpha_iN < alphas[iU]) {
					alphas[iU] = alpha_iN;
				}

			}
			else if ((uiN[iN][iU] - u0[iU]) < 0) {

				double alpha_iN = min(1.0, (uMIN[iU] - u0[iU]) / (uiN[iN][iU] - u0[iU]));

				if (alpha_iN < alphas[iU]) {
					alphas[iU] = alpha_iN;
				}

			}
			else {

				double alpha_iN = 1.0;

				if (alpha_iN < alphas[iU]) {
					alphas[iU] = alpha_iN;
				}
			}

		} // end loop for each state variable
	} // end loop for each node

		//out alpha
		// cout << "Alpha" << endl;
		// for (int i = 0; i < alphas.size(); i++) {
		//     cout <<std::setprecision(15)<< alphas[i] << " ";
		// cout << endl;
		// }


	// Scale the limeters
	for (int iU = 0; iU < 4; iU++) {

		L[iU] *= alphas[iU];

	}

	return L;

}


// vector<int> findAdjElem(vector<vector<double>> const& elem, vector<vector<double>> const& I2E, vector<vector<double>> const& B2E, vector<vector<int>> const& elemBounds, vector<vector<double>> const& bounds, vector<vector<double>> const& interiorFaces, vector<vector<int>> const& globalEdge, int iFaceGlobal) {

// 	// COMPUTES THE L and R elems for an input face (if no adj element exists due to being on a boundary, negative value is used)
// 	iGlobal2Local iLocal = iG2L(iFaceGlobal, bounds, interiorFaces);
// 	vector<int> elemLR = { -1, -1 };

// 	if (iLocal.isBound == false) { // interior edge
// 		elemLR[0] = int(I2E[iLocal.index][0] - 1); // Left Element
// 		elemLR[1] = int(I2E[iLocal.index][2] - 1); // Right Element
// 	}

// 	else { // boundary edge
// 		elemLR[0] = int(B2E[iLocal.index][0] - 1); // Left Element
// 		// No right element exists
// 	}

// 	return elemLR; // RETURNING L and R ELEMENT TO EDGE IN BASE-0 INDEXING

// }


// vector<int> findNeighbors(vector<vector<double>> const& elem, vector<vector<double>> const& I2E, vector<vector<double>> const& B2E, vector<vector<int>> const& elemBounds, vector<vector<double>> const& bounds, vector<vector<double>> const& interiorFaces, vector<vector<int>> const& globalEdge, int iElem) {
// 	// iElem is 0-based
// 	unordered_set<int> neighbors;

// 	// finding all edges surrounding an element
// 	vector<int> edges = { elemBounds[iElem][0], elemBounds[iElem][1], elemBounds[iElem][2] };

// 	// populating the neighbor set
// 	for (int iEdge = 0; iEdge < 3; iEdge++) {
// 		neighbors.insert(edges[iEdge]);
// 	}

// 	// Populating vector containing all neighbor indices
// 	vector<int> iNeighbors;
// 	iNeighbors.insert(iNeighbors.end(), neighbors.begin(), neighbors.end());

// 	// If we are on a boundary, some edges do not have an adjacent element, therefore indicate this with negative values
// 	if (iNeighbors.size() == 1) {
// 		iNeighbors.push_back(-1);
// 		iNeighbors.push_back(-1);
// 	}
// 	else if (iNeighbors.size() == 2) {
// 		iNeighbors.push_back(-1);
// 	}

// 	return iNeighbors;

// }

double computeEdgeLength(int iFace, vector<vector<double>> const& face, vector<vector<double>> const& nodes) {
	// Computes the length of an edge given the global edge index and the node coordinates
	double n1x = nodes[int(face[iFace][0]) - 1][0];
	double n1y = nodes[int(face[iFace][0]) - 1][1];
	double n2x = nodes[int(face[iFace][1]) - 1][0];
	double n2y = nodes[int(face[iFace][1]) - 1][1];

	double deltaL = sqrt(pow(n2x - n1x, 2) + pow(n2y - n1y, 2));
	return deltaL;
}


// Centroid of triangular element iElem.
inline Vector2d elementCentroid(vector<vector<double>> const& nodes, vector<vector<double>> const& elem, int iElem) {
	double x = (nodes[int(elem[iElem][0]) - 1][0] + nodes[int(elem[iElem][1]) - 1][0] + nodes[int(elem[iElem][2]) - 1][0]) / 3.0;
	double y = (nodes[int(elem[iElem][0]) - 1][1] + nodes[int(elem[iElem][1]) - 1][1] + nodes[int(elem[iElem][2]) - 1][1]) / 3.0;
	return { x, y };
}

// Everything both passes below need about one face. Global face indices run over
// the interior edges first and the boundary edges second, matching the numbering
// genE2F assigns - iG2L in processMesh.h decodes the same split.
struct faceGeom {
	int iElemL;             // left element, 0-based
	int iElemR;             // right element, 0-based; -1 on a boundary face
	bool isBound;
	int iFaceLocal;         // index into bounds/B2E/Bn, or interiorFaces/I2E/In
	vector<double> n;       // unit normal pointing OUT of iElemL
	Vector2d midpoint;
	double length;
};

// n is guaranteed to point out of iElemL: genIn orients In from the left element
// toward the right, and genBn orients Bn out of the domain. Both directions are
// verified by the per-cell closure sum(n_out * dl) = 0 that verification() checks,
// so callers can use f.n without re-deriving its sign.
inline faceGeom computeFaceGeom(meshData const& mesh, int iFaceGlobal) {
	faceGeom f;
	f.isBound    = iFaceGlobal >= mesh.niedge;
	f.iFaceLocal = f.isBound ? iFaceGlobal - mesh.niedge : iFaceGlobal;

	// reference, not a copy: this is per-face, and these tables are mesh-sized
	vector<vector<double>> const& face = f.isBound ? mesh.bounds : mesh.interiorFaces;

	if (f.isBound) {
		f.iElemL = int(mesh.B2E[f.iFaceLocal][0]) - 1;
		f.iElemR = -1;
		f.n = mesh.Bn[f.iFaceLocal];
	} else {
		f.iElemL = int(mesh.I2E[f.iFaceLocal][0]) - 1;
		f.iElemR = int(mesh.I2E[f.iFaceLocal][2]) - 1;
		f.n = mesh.In[f.iFaceLocal];
	}

	double n1x = mesh.nodes[int(face[f.iFaceLocal][0]) - 1][0];
	double n1y = mesh.nodes[int(face[f.iFaceLocal][0]) - 1][1];
	double n2x = mesh.nodes[int(face[f.iFaceLocal][1]) - 1][0];
	double n2y = mesh.nodes[int(face[f.iFaceLocal][1]) - 1][1];

	f.midpoint = { 0.5 * (n1x + n2x), 0.5 * (n1y + n2y) };
	f.length   = computeEdgeLength(f.iFaceLocal, face, mesh.nodes);

	return f;
}

// Riemann solver selection. The numbering matches the prompt FVM1st shows the user.
enum fluxOption { FLUX_ROE = 1, FLUX_RUSANOV = 2, FLUX_HLLE = 3 };

inline structFlux computeFlux(int opt, vector<double> const& uL, vector<double> const& uR, vector<double> const& n) {
	switch (opt) {
		case FLUX_ROE:     return roe(uL, uR, 1.4, n);
		case FLUX_RUSANOV: return rusanov(uL, uR, 1.4, n);
		case FLUX_HLLE:    return HLLE(uL, uR, 1.4, n);
	}
	// an unrecognised option would otherwise leave the flux struct uninitialised
	// and quietly poison every residual
	cerr << "ERROR: unknown flux option " << opt << " (expected 1 = Roe, 2 = Rusanov, 3 = HLLE)\n";
	exit(EXIT_FAILURE);
}

/* ===================================================================
 * Fast path: precomputed geometry + flat state arrays
 *
 * The mesh does not move, so everything the residual needs about a face is
 * computed once here instead of once per face per pass per RK stage: the
 * normal, the edge length, and - since the centroids are fixed too - the
 * reconstruction offset (face midpoint minus cell centroid) directly, so the
 * inner loop never touches mesh.nodes or mesh.elem at all.
 *
 * State and residual are flat arrays (stride 4 and 5) rather than
 * vector<vector<double>>: one contiguous block instead of nelem separate heap
 * allocations chased through a pointer table.
 * =================================================================== */

struct faceCache {
	int iElemL;
	int iElemR;        // -1 on a boundary face
	int isWall;        // only meaningful when iElemR < 0
	double nx, ny;     // unit normal, pointing out of iElemL
	double len;
	double dxL, dyL;   // face midpoint minus centroid of the left element
	double dxR, dyR;   // ... and of the right element (unused on a boundary)
};

struct meshGeom {
	vector<faceCache> faces;
	double uinf[4];    // freestream state, so the farfield ghost costs no trig
};

// Build once, right after buildMeshTopology, and pass to every residual call.
inline meshGeom buildMeshGeom(meshData const& mesh, double Minf, double alphaDeg) {
	meshGeom g;
	int nfaces = mesh.niedge + mesh.nbedge;
	g.faces.resize(nfaces);

	// Centroids are needed per face but there are only nelem of them, so they
	// are gathered once here and folded into the per-face offsets below.
	vector<Vector2d> cent(mesh.nelem);
	for (int e = 0; e < mesh.nelem; e++) cent[e] = elementCentroid(mesh.nodes, mesh.elem, e);

	for (int i = 0; i < nfaces; i++) {
		faceGeom f = computeFaceGeom(mesh, i);
		faceCache& c = g.faces[i];
		c.iElemL = f.iElemL;
		c.iElemR = f.isBound ? -1 : f.iElemR;
		c.isWall = (f.isBound && mesh.bounds[f.iFaceLocal][3] == 1) ? 1 : 0;
		c.nx = f.n[0];
		c.ny = f.n[1];
		c.len = f.length;
		c.dxL = f.midpoint[0] - cent[f.iElemL][0];
		c.dyL = f.midpoint[1] - cent[f.iElemL][1];
		if (c.iElemR >= 0) {
			c.dxR = f.midpoint[0] - cent[f.iElemR][0];
			c.dyR = f.midpoint[1] - cent[f.iElemR][1];
		} else {
			c.dxR = c.dyR = 0.0;
		}
	}

	vector<double> uinf = computeFreestreamState(Minf, alphaDeg);
	for (int k = 0; k < 4; k++) g.uinf[k] = uinf[k];
	return g;
}

// Scratch space for the gradient, owned by the caller so it is allocated once
// for the whole solve rather than twice per time step.
struct fvWorkspace {
	vector<double> grad;   // [nelem][4][x,y], flattened - stride 8 per element
	explicit fvWorkspace(int nelem) : grad(size_t(nelem) * 8, 0.0) {}
};

// Flat-array counterpart of computeBoundaryState: same wall projection and same
// farfield state, writing four doubles instead of returning a vector.
inline void boundaryStateFlat(const double* U_cell, bool isWall, const double* uinf,
							  double nx, double ny, double* Ub) {
	if (!isWall) { // farfield
		Ub[0] = uinf[0]; Ub[1] = uinf[1]; Ub[2] = uinf[2]; Ub[3] = uinf[3];
		return;
	}
	double rho = U_cell[0];
	double vx = U_cell[1] / rho, vy = U_cell[2] / rho;
	double vn = vx * nx + vy * ny;
	Ub[0] = rho;
	Ub[1] = rho * (vx - vn * nx);
	Ub[2] = rho * (vy - vn * ny);
	Ub[3] = U_cell[3];
}

// The flux is a template parameter, not a switch inside the loop: a function
// pointer known at compile time is a direct, inlinable call, so the flux body
// is folded into the face loop instead of costing a call per face.
template <double (*FluxCore)(const double*, const double*, double, double, double, double*)>
static void secondOrderResidualImpl(meshData const& mesh, meshGeom const& g,
									fvWorkspace& ws, const double* u, double* R) {
	const int nelem = mesh.nelem;
	double* grad = ws.grad.data();
	fill(ws.grad.begin(), ws.grad.end(), 0.0);
	fill(R, R + size_t(nelem) * 5, 0.0);

	// Pass 1: grad(u)_i = (1/A_i) * sum over faces of uhat * n_out * dl
	for (size_t i = 0; i < g.faces.size(); i++) {
		faceCache const& f = g.faces[i];
		const double* uL = u + size_t(f.iElemL) * 4;
		double ghost[4];
		const double* uR;
		if (f.iElemR < 0) {
			boundaryStateFlat(uL, f.isWall != 0, g.uinf, f.nx, f.ny, ghost);
			uR = ghost;
		} else {
			uR = u + size_t(f.iElemR) * 4;
		}

		double* gL = grad + size_t(f.iElemL) * 8;
		if (f.iElemR >= 0) {
			double* gR = grad + size_t(f.iElemR) * 8;
			for (int k = 0; k < 4; k++) {
				double c = 0.5 * (uL[k] + uR[k]) * f.len;
				double cx = c * f.nx, cy = c * f.ny;
				gL[k * 2] += cx; gL[k * 2 + 1] += cy;
				// the same normal points INTO the right element, hence the sign flip
				gR[k * 2] -= cx; gR[k * 2 + 1] -= cy;
			}
		} else {
			for (int k = 0; k < 4; k++) {
				double c = 0.5 * (uL[k] + uR[k]) * f.len;
				gL[k * 2] += c * f.nx; gL[k * 2 + 1] += c * f.ny;
			}
		}
	}

	// Dividing each gradient by its element area
	for (int e = 0; e < nelem; e++) {
		double* ge = grad + size_t(e) * 8;
		double A = mesh.area[e];
		for (int k = 0; k < 8; k++) ge[k] /= A;
	}

	// Pass 2: reconstruct to each face midpoint and accumulate the flux residual
	for (size_t i = 0; i < g.faces.size(); i++) {
		faceCache const& f = g.faces[i];
		const double* gL = grad + size_t(f.iElemL) * 8;
		const double* cL = u + size_t(f.iElemL) * 4;
		double uLr[4], uRr[4];
		for (int k = 0; k < 4; k++)
			uLr[k] = cL[k] + gL[k * 2] * f.dxL + gL[k * 2 + 1] * f.dyL;

		if (f.iElemR >= 0) {
			const double* gR = grad + size_t(f.iElemR) * 8;
			const double* cR = u + size_t(f.iElemR) * 4;
			for (int k = 0; k < 4; k++)
				uRr[k] = cR[k] + gR[k * 2] * f.dxR + gR[k * 2 + 1] * f.dyR;
		} else {
			// Built from the already-reconstructed uLr, so the wall tangency holds
			// at the point the flux actually sees.
			boundaryStateFlat(uLr, f.isWall != 0, g.uinf, f.nx, f.ny, uRr);
		}

		double F[4];
		double s = FluxCore(uLr, uRr, 1.4, f.nx, f.ny, F);

		double* RL = R + size_t(f.iElemL) * 5;
		for (int k = 0; k < 4; k++) RL[k] += F[k] * f.len;
		RL[4] += s * f.len;

		if (f.iElemR >= 0) {
			double* RR = R + size_t(f.iElemR) * 5;
			for (int k = 0; k < 4; k++) RR[k] -= F[k] * f.len;
			RR[4] += s * f.len;
		}
	}
}

// R must have room for nelem*5 doubles, u must hold nelem*4.
inline void secondOrderResidual(meshData const& mesh, meshGeom const& g, fvWorkspace& ws,
								int opt, const double* u, double* R, string const& limiterType) {
	if (limiterType != "NONE") {
		cerr << "ERROR: limiterType \"" << limiterType << "\" is not implemented yet.\n"
			 << "       Only \"NONE\" (unlimited reconstruction) is available; the BJ and LCD\n"
			 << "       limiters exist as functions but are not called from this driver.\n";
		exit(EXIT_FAILURE);
	}
	switch (opt) {
		case FLUX_ROE:     secondOrderResidualImpl<roeCore>(mesh, g, ws, u, R);     return;
		case FLUX_RUSANOV: secondOrderResidualImpl<rusanovCore>(mesh, g, ws, u, R); return;
		case FLUX_HLLE:    secondOrderResidualImpl<hlleCore>(mesh, g, ws, u, R);    return;
	}
	cerr << "ERROR: unknown flux option " << opt << " (expected 1 = Roe, 2 = Rusanov, 3 = HLLE)\n";
	exit(EXIT_FAILURE);
}

// L1 norm over the four residual components of a flat [nelem][5] residual.
inline double computeL1ResidualNorm(const double* residuals, int nelem) {
	double L1ResidualNorm = 0;
	for (int iElem = 0; iElem < nelem; iElem++) {
		for (int iU = 0; iU < 4; iU++) {
			L1ResidualNorm += abs(residuals[size_t(iElem) * 5 + iU]);
		}
	}
	return L1ResidualNorm;
}

// 2nd Order Finite Volume Driver - REFERENCE IMPLEMENTATION
//
// Superseded by secondOrderResidual above, which is what rk2 calls. This one is
// kept because it is the readable statement of the scheme and the baseline the
// fast path is checked against (see profile_fast.cpp); it is not on the hot path.
//
// Two sweeps over the faces: the first accumulates the Green-Gauss gradient in
// every cell, the second reconstructs the state to each face midpoint and
// accumulates the flux residual. Returns [nelem][5] - the four residual
// components, then the edge-length-weighted wave speed tally that rk2 divides
// into the cell area for its local time step.
vector<vector<double>> secondOrderFV(meshData const& mesh, int opt, vector<vector<double>>& u,
									 double const& Minf, double alphaDeg, string limiterType) {

	// barthJespersen and computeL_LCD above are not wired into this driver, so a
	// "BJ" or "MP" run would silently produce an unlimited reconstruction. Refuse
	// it rather than let unlimited results be written up as limited ones.
	if (limiterType != "NONE") {
		cerr << "ERROR: limiterType \"" << limiterType << "\" is not implemented in secondOrderFV yet.\n"
		     << "       Only \"NONE\" (unlimited reconstruction) is available; the BJ and LCD\n"
		     << "       limiters exist as functions but are not called from this driver.\n";
		exit(EXIT_FAILURE);
	}

	int nelem = mesh.nelem;
	int nfaces = mesh.niedge + mesh.nbedge;

	// Pass 1: grad(u)_i = (1/A_i) * sum over faces of uhat * n_out * dl
	vector<vector<Vector2d>> grad_u(nelem, vector<Vector2d>(4, Vector2d::Zero()));

	for (int iFace = 0; iFace < nfaces; iFace++) {
		faceGeom f = computeFaceGeom(mesh, iFace);
		Vector2d nOutL(f.n[0], f.n[1]);

		vector<double> const& uL = u[f.iElemL];
		vector<double> uR = f.isBound
			? computeBoundaryState(uL, mesh.bounds[f.iFaceLocal][3] == 1, Minf, alphaDeg, mesh.Bn, f.iFaceLocal)
			: u[f.iElemR];

		for (int iU = 0; iU < 4; iU++) {
			Vector2d contribution = 0.5 * (uL[iU] + uR[iU]) * f.length * nOutL;
			grad_u[f.iElemL][iU] += contribution;
			// the same normal points INTO the right element, hence the sign flip
			if (!f.isBound) grad_u[f.iElemR][iU] -= contribution;
		}
	}

	// Dividing each grad_u by its element area
	for (int iElem = 0; iElem < nelem; iElem++) {
		for (int iU = 0; iU < 4; iU++) {
			grad_u[iElem][iU] /= mesh.area[iElem];
		} // end for iU
	} // end for iElem

	// Initialize Residual and Wave Speed on Each Cell to be Zero
	vector<vector<double>> residuals(nelem, vector<double>(5, 0)); // residual for each state and then the wave speed

	for (int iFace = 0; iFace < nfaces; iFace++) {
		faceGeom f = computeFaceGeom(mesh, iFace);
		Vector2d centroidL = elementCentroid(mesh.nodes, mesh.elem, f.iElemL);

		// u_L extrapolated from the left cell centroid out to the face midpoint
		vector<double> uL = u[f.iElemL];
		for (int iU = 0; iU < 4; iU++) {
			uL[iU] += grad_u[f.iElemL][iU].dot(f.midpoint - centroidL);
		}

		vector<double> uR;
		if (!f.isBound) {
			Vector2d centroidR = elementCentroid(mesh.nodes, mesh.elem, f.iElemR);
			uR = u[f.iElemR];
			for (int iU = 0; iU < 4; iU++) {
				uR[iU] += grad_u[f.iElemR][iU].dot(f.midpoint - centroidR);
			}
		} else {
			// Built from the already-reconstructed uL, so the wall tangency holds at
			// the point the flux actually sees. Extrapolating a second time out to a
			// mirrored ghost centroid would undo it - that offset is exactly
			// -(midpoint - centroidL).
			uR = computeBoundaryState(uL, mesh.bounds[f.iFaceLocal][3] == 1, Minf, alphaDeg, mesh.Bn, f.iFaceLocal);
		}

		structFlux output = computeFlux(opt, uL, uR, f.n);

		// R_i = sum over faces of Fhat . n_out * dl. f.n points out of L, so this
		// face adds to L and subtracts from R; a boundary face has no right element
		// and only adds. The wave speed is a magnitude, so it accumulates on both.
		for (int j = 0; j < 4; j++) {
			residuals[f.iElemL][j] += output.F[j] * f.length;
		}
		residuals[f.iElemL][4] += output.s_mag * f.length;

		if (!f.isBound) {
			for (int j = 0; j < 4; j++) {
				residuals[f.iElemR][j] -= output.F[j] * f.length;
			}
			residuals[f.iElemR][4] += output.s_mag * f.length;
		}

	} // end for iFace

	return residuals;

}

#endif /* solver_h */
