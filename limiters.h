#include <vector>
#include <algorithm>
#include <iostream>


// Obtaining all elements neighboring elemL
vector<int> iNeighbor;
vector<int> E2F_i = mesh.E2F[iElemL]; // global face indices for the current element
		
    // For each face, find the two elements bounding the face, 
    // and push the element which is not elemL into iNeighbor
for (int iF = 0; iF < 3; iF++) {
    vector<int> adjElems = findAdjElem(mesh.elem, mesh.I2E, mesh.B2E, mesh.elemBounds, mesh.bounds, mesh.interiorFaces, globalEdge, globalFaceIndices[iF]);
    auto it = find(adjElems.begin(), adjElems.end(), iElemL);

    adjElems.erase(it);
    if (adjElems.size() == 1 && adjElems[0] >= 0) {
        // !adjElems.empty()
        iNeighbor.push_back(adjElems[0]);

    }
    else {

        iGlobal2Local neighborBoundFace = iG2L(globalFaceIndices[iF], mesh.bounds, mesh.interiorFaces);
        int iNeighborFaceLocal = neighborBoundFace.index;
        if (mesh.bounds[iNeighborFaceLocal][3] == 0) { // farfield
            iNeighbor.push_back(-1);
        }
        else { // wall
            iNeighbor.push_back(-2);
        }


    }

		}

		// // Obtaining all edges bounding an element and local face numbering for these edges
		// vector<int> iFaces = { 0,1,2 }; // adjacent elements found using elemBounds, therefore the edges were indexed in order of the local face index

		// // Initializing the normal vector pointing from L to R
		// Vector2d n;

		if (isBoundFace == false) { // interior edge

			n = { mesh.In[iFaceLocal][0], mesh.In[iFaceLocal][1] };

		}
		else { // boundary edge

			n = { mesh.Bn[iFaceLocal][0], mesh.Bn[iFaceLocal][1] };

		}

		// Gradient Calculation
		vector<Vector2d> L;

		if (limiterType == "BJ") {

			L = barthJespersen(mesh.nodes, mesh.elem, mesh.area, U, iElemL, iNeighbor, Minf, alphaDeg, mesh.Bn, mesh.elemBounds, mesh.bounds, mesh.interiorFaces, iFaces);

		}
		else if (limiterType == "MP") {

			vector<Vector2d> L_input = computeL(mesh.nodes, mesh.elem, U, iElemL, iNeighbor, iFaces, Minf, alphaDeg, mesh.Bn, mesh.bounds, mesh.interiorFaces, mesh.elemBounds);
			L = computeL_LCD(L_input, U, mesh.nodes, mesh.elem, iElemL, iNeighbor, iFaces, Minf, alphaDeg, mesh.Bn, mesh.elemBounds, mesh.bounds, mesh.interiorFaces);

		}
		else if (limiterType == "NONE") {

			L = computeL(mesh.nodes, mesh.elem, U, iElemL, iNeighbor, iFaces, Minf, alphaDeg, mesh.Bn, mesh.bounds, mesh.interiorFaces, mesh.elemBounds);

		}
		else {

			cout << "Valid Limiter Type not Used\n";
			return vector<vector<double>>(0);

		}

		if (iElemL == 446 || iElemR == 446) {
			int stop = 0;
		}

		// // Compute face length
		double delta_l = computeEdgeLength(iFaceGlobal, globalEdge, mesh.nodes);

		// // Find u_hat (the average of the L and R cell averages)
		vector<double> UL_limiting = U[iElemL];

		vector<double> UR_limiting;
		UR_limiting.reserve(4);
		// If no UR exists, need to create a ghost state
		if (isBoundFace == true) {

			bool isWall = false;

			if (mesh.bounds[iFaceLocal][3] == 1) { // 1 = wall, 0 = farfield (from the .gri group title)

				isWall = true;

			}
			//                else{
			//                    int stop = 0;
			//
			//                }

			UR_limiting = computeBoundaryState(UL_limiting, isWall, Minf, alphaDeg, mesh.Bn, iFaceLocal);

		}
		else {
			UR_limiting = U[iElemR];
		}

		vector<double> Uhat;
		Uhat.reserve(UL_limiting.size());

		// // Add u_hat*n*delta_l to gradient value in vector for L elem and subtract this for the right element (if existing)

		for (int iU = 0; iU < UL_limiting.size(); iU++) {
			Uhat.emplace_back(0.5 * (UL_limiting[iU] + UR_limiting[iU]));
		}

		// // Adding u_hat*n*delta_l to elemL and subtracting from elemR
		for (int iU = 0; iU < UL_limiting.size(); iU++) {
			Vector2d gradU_add = { n[0] * Uhat[iU] * delta_l,n[1] * Uhat[iU] * delta_l };

			if (iElemR != -1) { // interior face
				grad_u[iElemL][iU] += gradU_add;
			}
			else { // boundary face, normal points out of element, therefore
				gradU_add = -gradU_add;
				grad_u[iElemL][iU] += gradU_add;
			}

			if (iLocal.isBound == false) { // if right element exists
				grad_u[iElemR][iU] -= gradU_add;
			}

		} // end for iU*/