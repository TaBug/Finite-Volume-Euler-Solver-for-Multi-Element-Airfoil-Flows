#ifndef processMesh_h
#define processMesh_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <numeric>
#include <cmath>
#include <filesystem>
#include <cstdlib>
using namespace std;

// Decides the boundary condition from a .gri boundary group title. The farfield
// box is named for its sides; every other group (the airfoil elements: main,
// slat, flap) is a solid wall. Unrecognised titles fall through to wall, which
// is the safer default and matches how the group test used to behave.
inline bool isWallTitle(const string &title) {
	return !(title == "bot" || title == "right" || title == "top" || title == "left" ||
	         title == "farfield" || title == "Farfield" || title == "freestream");
}

// STRUCTURE TO HOLD MESH DATA
struct meshData {
	int nnode = 0;
	int nbedge = 0;
	int nelem = 0;
	int niedge = 0; // interior edges, derived: (3*nelem - nbedge) / 2
	vector<vector<double>> nodes; // [nodeID][x, y]
	vector<vector<double>> elem; // [elemID][node1, node2, node3]
	vector<string> bTitles; // .gri boundary group titles, indexed by group - 1
	// [node1, node2, bGroup, isWall] - isWall is 1 for a solid wall, 0 for
	// farfield, resolved from the group title by isWallTitle() when the mesh is
	// read, so the solver never has to hard-code group numbers
	vector<vector<double>> bounds; // [boundID][node1, node2, bGroup, isWall]
	vector<vector<double>> interiorFaces; // [faceID][node1, node2, leftElem, rightElem] - left/right element
	vector<vector<double>> I2E; // [interiorFaceID][leftElem, leftFaceLocal, rightElem, rightFaceLocal]
	vector<vector<double>> B2E; // [boundaryFaceID][elem, faceLocal, bGroup]
	vector<vector<int>> E2F; // [elemID][face1, face2, face3] - global edge indices for each element
	vector<vector<double>> In; // [interiorFaceID][nx, ny] - normal vector for interior faces, pointing from left to right element
	vector<vector<double>> Bn;	// [boundaryFaceID][nx, ny] - normal vector for boundary faces
	vector<double> area;
};

// Returns the full path of the named mesh file, which lives in the "gri"
// subfolder of the working directory. Run the program from the project root.
inline string findGriFile(const string &fileName) {
	filesystem::path path = filesystem::current_path() / "gri" / fileName;

	// stop here rather than hand back a path that cannot be opened: the readers
	// downstream do not check, and would leave the solve running on empty data
	if (!filesystem::exists(path)) {
		cerr << "ERROR: mesh file not found: " << path.string() << "\n";
		exit(EXIT_FAILURE);
	}

	return path.string();
}

// READ GMSH .BDF OUTPUT FILE INTO NODE, BOUNDARY, AND ELEMENT MATRICES
// User input a Gmsh ".bdf" file type
/// PASS IN THE PARAMETERS WE NEED BY REFERENCE /////
void readGriFile (string fileName, struct meshData &mesh) {
    // Open the file for reading
    ifstream file(fileName);
	
	int discard; // Dummy variable for parameters we do not care about

    // If there is a file in the directory matching the input parameter "fileName", open it
    if (file.is_open()) {
        string currLine; // String to store each line of the file
		string type; // INITIALIZE DATA STRUCTURE SIZES

		// Read the first line of the file to get the number of nodes, elements, and dimensions
		getline(file, currLine);
		stringstream ss(currLine);
		ss >> mesh.nnode >> mesh.nelem >> discard;

		// Initialize the mesh node and element data structures
		mesh.nodes.resize(mesh.nnode, vector<double>(2));
		mesh.elem.resize(mesh.nelem, vector<double>(3));

		// Read the node data from the file
		for (int i = 0; i < mesh.nnode; i++) {
			getline(file, currLine);
			stringstream ss(currLine);
			ss >> mesh.nodes[i][0] >> mesh.nodes[i][1];
		}

		// Read the number of boundary groups from the file
		int nBGroup;
		getline(file, currLine);
		stringstream ssGroup(currLine);
		ssGroup >> nBGroup;

		// Read the boundary data from the file
		mesh.bounds.clear();
		mesh.bTitles.clear();
		for (int i = 0; i < nBGroup; i++) {
			int nBFace;
			string bTitle;
			getline(file, currLine);
			stringstream ss(currLine);
			// nBFace, nodes per face (2 = linear), then the group title. The
			// title is a name ("main", "bot"), not a number, so it has to be
			// read as a string - reading it as an int silently yields 0.
			ss >> nBFace >> discard >> bTitle;

			mesh.bTitles.push_back(bTitle);
			int bGroup = i + 1; // groups are numbered by their order in the file
			double isWall = isWallTitle(bTitle) ? 1.0 : 0.0;

			// Read the boundary face data for the current boundary group
			for (int j = 0; j < nBFace; j++) {
				getline(file, currLine);
				stringstream ss(currLine);
				int node1, node2;
				ss >> node1 >> node2;
				mesh.bounds.push_back({static_cast<double>(node1), static_cast<double>(node2),
									   static_cast<double>(bGroup), isWall});
			}
		}
		mesh.nbedge = int(mesh.bounds.size()); // only known once every group is read

		// Read the number of elements and discard the other two values. The
		// basis name ("TriLagrange") is text, so it is read as a string.
		string basis;
		getline(file, currLine);
		stringstream ssElem(currLine);
		ssElem >> discard >> discard >> basis;

		// Read the element data from the file
		for (int i = 0; i < mesh.nelem; i++) {
			getline(file, currLine);
			stringstream ss(currLine);
			int node1, node2, node3;
			ss >> node1 >> node2 >> node3;
			mesh.elem[i] = {static_cast<double>(node1), static_cast<double>(node2), static_cast<double>(node3)};
		};

        // Close the file
        file.close();
    } else {
        // stop rather than return a half-filled struct: nothing downstream
        // checks, so the solve would run on empty matrices
        cerr << "ERROR: unable to open mesh file " << fileName << "\n";
        exit(EXIT_FAILURE);
    }

	// Every triangle has 3 edges, each interior edge shared by 2 elements and
	// each boundary edge by 1, so 3*nelem - nbedge = 2*niedge must come out
	// positive and even. Catches a truncated or malformed file here rather than
	// as a bad allocation later.
	if (mesh.nnode == 0 || mesh.nelem == 0 ||
		3 * mesh.nelem - mesh.nbedge <= 0 || (3 * mesh.nelem - mesh.nbedge) % 2 != 0) {
		cerr << "ERROR: inconsistent mesh in " << fileName << " (nnode=" << mesh.nnode
		     << ", nbedge=" << mesh.nbedge << ", nelem=" << mesh.nelem << ")\n";
		exit(EXIT_FAILURE);
	}

	mesh.niedge = (3 * mesh.nelem - mesh.nbedge) / 2;

	// The derived connectivity is built by buildMeshTopology(), which is defined
	// below the generators it calls - call it after this function returns.
}

// CREATE A VECTOR WHICH CONTAINS THE AMOUNT OF EDGES IN EACH BOUNDARY GROUP
vector<int> boundSizes(vector<vector<double>> &bounds,int nbedge){
	
	unordered_set<int> unique_values;
	
	for(int i=0;i < nbedge; i++){
		unique_values.insert(bounds[i][2]);
	}
	
	size_t nBGroup = unique_values.size();
	vector<int> nBFace(nBGroup,0); // Vector which reports the # of bound faces in each bound group. Vector index represents bGroup
	
	for(int i=0;i < nbedge; i++){
		nBFace[bounds[i][2] - 1]++; // incrementing group index by one
	}
	
	return nBFace;
}


// TASK 2
/* Construct a matrix that maps the global face index to the indices of the end nodes, left element index, and right element index
	INPUTS: niedge = total number of edges
			nelem = total number of elements
			bounds = matrix that stores the boudary edges
			elem = matrix that stores the vertex indices of the elements 
	OUTPUTS: interiorFaces = matrix that stores the indices of two end nodes, left element index, and right element index (# of faces x 4) */
vector<vector<double>> genInteriorFaceVec(struct meshData &mesh) {
	vector<vector<double>> interiorFaces;
	int facePairs[3][2] = { {0, 1}, {1, 2}, {2, 0} }; // Pairs of nodes that make up each face of a triangle
	
	// Creating a set of all nodes that lie on a boundary
	unordered_set<double> bNodeSet;
	for(int i = 0; i < mesh.nbedge; i++) {
		bNodeSet.insert(mesh.bounds[i][0]);
		bNodeSet.insert(mesh.bounds[i][1]);
	}
	
	// Parsing through all faces and storing those on the interior into the interiorFaces vector
	for(int iElem = 0; iElem < mesh.nelem; iElem++) {
		for (int iFace = 0; iFace < 3; iFace++) {
			// Parsing through each face of the current element
			double faceNode1 = mesh.elem[iElem][facePairs[iFace][0]];
			double faceNode2 = mesh.elem[iElem][facePairs[iFace][1]];

			// Check to see if both nodes do not lie on a boundary
			if (bNodeSet.find(faceNode1) == bNodeSet.end() || 
				bNodeSet.find(faceNode2) == bNodeSet.end()) {
				// Ensure face is not already stored
				bool inVec = false; // true if face is stored, otherwise false
				int indexMatchingFace = -1; // face global index


				for(int i = 0; i < interiorFaces.size(); i++) {
					// if the face is already stored, record the global index
					if( (faceNode1 == interiorFaces[i][0] && faceNode2 == interiorFaces[i][1]) || 
						(faceNode2 == interiorFaces[i][0] && faceNode1 == interiorFaces[i][1]) ) {
						indexMatchingFace = i;
						inVec = true;
					}
				}

				// If boundary face not yet stored, insert it into interiorFaces vector
				if (inVec == false) {
					interiorFaces.push_back({faceNode1, faceNode2, double(iElem) + 1});
				} else { // If boundary face is already stored, add the second element that is adjacent to face
					interiorFaces[indexMatchingFace].push_back(iElem + 1);
				}

				// Reset the flag
				inVec = false;
			}
		}
	}
	return interiorFaces;
}

vector<vector<double>> genI2E (meshData &mesh) {
	
	vector<vector<double>> I2E(mesh.niedge, vector<double>(4));
	
	int intEdgeNum = int(mesh.interiorFaces.size());
	
	for(int i = 0; i < intEdgeNum; i++){
		double node1 = mesh.interiorFaces[i][0];
		double node2 = mesh.interiorFaces[i][1];
		double elemL = mesh.interiorFaces[i][2];
		double elemR = mesh.interiorFaces[i][3];
		// Calculate local face number for elemL
		double faceL = -1;
		vector<double> elemNodesL({mesh.elem[elemL-1][0],mesh.elem[elemL-1][1],mesh.elem[elemL-1][2]});
		for(int j = 0; j < 3; j++){
			if(elemNodesL[j] != node1 && elemNodesL[j] != node2){
				// faceL = elemNodesL[j];
				faceL = j+1;
				break;
			}
		}
		// Calculate local face number for elemR
		double faceR = -1;
		vector<double> elemNodesR({mesh.elem[elemR-1][0],mesh.elem[elemR-1][1],mesh.elem[elemR-1][2]});
		for(int j = 0; j < 3; j++){
			if(elemNodesR[j] != node1 && elemNodesR[j] != node2){
				// faceR = elemNodesR[j];
				faceR = j+1;
				break;
			}
		}
		
		I2E[i][0] = elemL;
		I2E[i][1] = faceL;
		I2E[i][2] = elemR;
		I2E[i][3] = faceR;
		
	}
	
	return I2E;
}

// NEEDS TO VERIFY
vector<vector<double>> genB2E(meshData &mesh){
	vector<vector<double>> B2E(mesh.nbedge, vector<double>(4));

	for(int i = 0; i < mesh.nbedge; i++){
		double n1 = mesh.bounds[i][0]; // Nodes on boundary
		double n2 = mesh.bounds[i][1]; // Nodes on boundary
		// Saving bgroup
		B2E[i][2] = mesh.bounds[i][2];
		// Saving the wall/farfield flag resolved from the group title
		B2E[i][3] = mesh.bounds[i][3];
		// find element containing nodes
		for(int j = 0; j < mesh.nelem; j++){
			
			double val1 = mesh.elem[j][0];
			double val2 = mesh.elem[j][1];
			double val3 = mesh.elem[j][2];
			
			if ((val1 == n1 && val2 == n2) || (val1 == n2 && val2 == n1) 
				|| (val2 == n1 && val3 == n2) || (val2 == n2 && val3 == n1) 
				|| (val1 == n1 && val3 == n2) || (val1 == n2 && val3 == n1)) {
				
				B2E[i][0] = double(j + 1);
				
				if(val1 != n1 && val1 != n2){
					// B2E[i][1] = elem[j][0];
					B2E[i][1] = 1;
				}
				else if(val2 != n1 && val2 != n2){
					// B2E[i][1] = elem[j][1];
					B2E[i][1] = 2;
				}
				else{
					// B2E[i][1] = elem[j][2];
					B2E[i][1] = 3;
				}
				
				break;
				
			}
		}
	}
	
	return B2E;
}

vector<vector<int>> genE2F(meshData const &mesh) {
    vector<vector<int>> E2F(mesh.nelem, vector<int>(3, 0));
    vector<vector<double>> I2E = mesh.I2E;
    vector<vector<double>> B2E = mesh.B2E;

    // Start with interior edges
    for(int i = 0; i < mesh.niedge; i++) {
        E2F[int(I2E[i][0]) - 1][int(I2E[i][1]) - 1] = i;
        E2F[int(I2E[i][2]) - 1][int(I2E[i][3])-1] = i;  
    }
    
    // Now do boundary edges
    for(int i = 0; i < mesh.nbedge; i++){
        E2F[int(B2E[i][0])-1][int(B2E[i][1])-1] = i + mesh.niedge;
    }

    return E2F;
}

// Result of splitting a global edge index back into the table it came from.
struct iGlobal2Local {
    bool isBound;   // true  -> index into bounds / B2E / Bn
                    // false -> index into interiorFaces / I2E / In
    int index;
};

// Inverse of the numbering genE2F assigns above: interior edges occupy
// [0, niedge) and boundary edges occupy [niedge, niedge + nbedge). Keep this
// next to genE2F - it is the only decoder of that layout, and the two silently
// disagree if either is renumbered on its own.
inline iGlobal2Local iG2L(int iGlobal, vector<vector<double>> const &bounds, vector<vector<double>> const &interiorFaces){

    int niedge = int(interiorFaces.size());

    if(iGlobal < niedge){
        return { false, iGlobal };
    }
    return { true, iGlobal - niedge };

}

vector<vector<double>> genIn(int niedge, vector<vector<double>> &interiorFaces, vector<vector<double>> &elem, vector<vector<double>> &nodes, vector<vector<double>> &I2E){
	
	vector<vector<double>> In(niedge, vector<double>(2));
	
	for(int i = 0; i < niedge; i++){
	
		double n1 = interiorFaces[i][0], n2 = interiorFaces[i][1]; // Nodes on interior face

		// The vertex of the LEFT element that is not on this face. Local face k of
		// a triangle is opposite local node k, so I2E's leftFaceLocal indexes it
		// directly. Pointing the normal away from this vertex is what makes In run
		// from the left element to the right one, as meshData documents.
		size_t iFaceLocal = size_t(I2E[i][1]) - 1;          // leftFaceLocal
		size_t iElem      = size_t(interiorFaces[i][2]);    // leftElem
		double fLocal = elem[iElem - 1][iFaceLocal];

		size_t ifLocal = fLocal - 1; // index of local node
		// Obtain coordinate of local node and boundary node
		double x1 = nodes[n1 - 1][0], x2 = nodes[n2 - 1][0];
		double y1 = nodes[n1 - 1][1], y2 = nodes[n2 - 1][1];
		double l = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
		double xLocal = nodes[ifLocal][0], yLocal = nodes[ifLocal][1];

		// Unit tangent along the face, n1 -> n2
		double xComponent = (1/l) * (x2 - x1);
		double yComponent = (1/l) * (y2 - y1);

		// Rotate the tangent 90 degrees to get a unit normal, then flip it if it
		// points back toward the left element's opposite vertex. Doing the
		// orientation with a dot product rather than a slope comparison keeps this
		// correct for vertical and horizontal faces, with no division by (x2 - x1).
		double xNorm =  yComponent;
		double yNorm = -xComponent;

		if (xNorm * (xLocal - x1) + yNorm * (yLocal - y1) > 0.0) {
			xNorm = -xNorm;
			yNorm = -yNorm;
		}

		In[i] = vector<double> {xNorm, yNorm};

	}
	
	return In;
}

vector<vector<double>> genBn(int nbedge, vector<vector<double>> &bounds, vector<vector<double>> &nodes, vector<vector<double>> &elem, vector<vector<double>> &B2E){
	
	vector<vector<double>> Bn(nbedge, vector<double>(2));
	
	for(int i = 0; i < nbedge; i++){
		// Nodes on  boundary
		double n1 = bounds[i][0], n2 = bounds[i][1];

		// The vertex of the boundary element that is not on this face. Local face k
		// of a triangle is opposite local node k, so B2E's faceLocal indexes it.
		size_t iFaceLocal = size_t(B2E[i][1]) - 1;
		size_t iElem = size_t(B2E[i][0]);
		double fLocal = elem[iElem - 1][iFaceLocal];
		size_t ifLocal = fLocal - 1; // index of local node
		// Obtain coordinate of local node and boundary node
		double x1 = nodes[n1 - 1][0], x2 = nodes[n2 - 1][0];
		double y1 = nodes[n1 - 1][1], y2 = nodes[n2 - 1][1];
		double xLocal = nodes[ifLocal][0], yLocal = nodes[ifLocal][1];

		double magnitude = sqrt(pow(x2-x1, 2) + pow(y2-y1,2));

		// Unit tangent along the face, n1 -> n2
		double xComponent = (1/magnitude)*(x2-x1);
		double yComponent = (1/magnitude)*(y2-y1);

		// Rotate the tangent 90 degrees, then flip it if it points back toward the
		// element's opposite vertex, so Bn ends up pointing OUT of the domain. The
		// slope-comparison version this replaces produced inward normals on every
		// boundary face, which the solver was compensating for by negating the
		// boundary contribution at the call site.
		double xNorm =  yComponent;
		double yNorm = -xComponent;

		if (xNorm * (xLocal - x1) + yNorm * (yLocal - y1) > 0.0) {
			xNorm = -xNorm;
			yNorm = -yNorm;
		}

		Bn[i] = vector<double> {xNorm,yNorm};

	}
	
	return Bn;
}

vector<double> genArea(int nelem, vector<vector<double>> &elem, vector<vector<double>> &nodes){
	vector<double> AREA(nelem);
	
	for(int i = 0; i < nelem; i++){
		
		double n1, n2, n3;
		double x1, y1, x2, y2, x3, y3;
		double a, b, c, s;
		
		n1 = elem[i][0];
		n2 = elem[i][1];
		n3 = elem[i][2];
		
		x1 = nodes[n1 - 1][0];
		x2 = nodes[n2 - 1][0];
		x3 = nodes[n3 - 1][0];
		
		y1 = nodes[n1 - 1][1];
		y2 = nodes[n2 - 1][1];
		y3 = nodes[n3 - 1][1];
		
		// Calculate the lengths of the three sides
		a = sqrt(pow(x2 - x3, 2) + pow(y2 - y3, 2));
		b = sqrt(pow(x1 - x3, 2) + pow(y1 - y3, 2));
		c = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
		
		// Calculate the semi-perimeter of the triangle
		s = (a + b + c) / 2;
		
		AREA[i] = sqrt(s * (s - a) * (s - b) * (s - c));;

	}

	return AREA;
}

// meshData overloads for the three generators that still take explicit
// arguments, so every generator can be called uniformly as gen*(mesh)
inline vector<vector<double>> genIn(meshData &mesh) {
	return genIn(mesh.niedge, mesh.interiorFaces, mesh.elem, mesh.nodes, mesh.I2E);
}

inline vector<vector<double>> genBn(meshData &mesh) {
	return genBn(mesh.nbedge, mesh.bounds, mesh.nodes, mesh.elem, mesh.B2E);
}

inline vector<double> genArea(meshData &mesh) {
	return genArea(mesh.nelem, mesh.elem, mesh.nodes);
}

// Builds every derived connectivity structure from the raw mesh read by
// readGriFile. Defined here because it calls the generators above; the order of
// the calls matters, since I2E feeds In and B2E feeds Bn.
inline void buildMeshTopology(meshData &mesh) {
	mesh.interiorFaces = genInteriorFaceVec(mesh);
	mesh.I2E = genI2E(mesh);
	mesh.B2E = genB2E(mesh);
	mesh.E2F = genE2F(mesh);
	mesh.In = genIn(mesh);
	mesh.Bn = genBn(mesh);
	mesh.area = genArea(mesh);
}

// VERIFICATION
// Calculates the error on each element for mesh verification
vector<double> verification(int nelem, vector<vector<double>> &I2E, vector<vector<double>> &B2E, vector<vector<double>> &In, vector<vector<double>> &Bn, vector<vector<double>> &nodes, vector<vector<double>> &elem, vector<vector<double>> &bounds){
	
	vector<vector<double>> error(nelem, vector<double>(2));
	vector<double> err_mag(nelem);
	int I2E_size = int(I2E.size());
	int B2E_size = int(B2E.size());

	// loop through interior edges
	for (int i = 0; i < I2E_size; i++){
		int elemL = I2E[i][0] - 1;
		int elemR = I2E[i][2] - 1;
		int face = I2E[i][1];
		double length = -1;
		vector<double> coord1(2);
		vector<double> coord2(2);
		

		// calculate length of face
		if (face == 1){
			int node1 = elem[elemL][1] - 1;
			int node2 = elem[elemL][2] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		else if (face == 2){
			int node1 = elem[elemL][0] - 1;
			int node2 = elem[elemL][2] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		else if (face == 3){
			int node1 = elem[elemL][0] - 1;
			int node2 = elem[elemL][1] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		
		// Compute error of interior edges with length and normal vectors. In now
		// points out of elemL (left -> right), so it is elemL that accumulates +In
		// and elemR that accumulates -In; the closure sum over each element is then
		// sum(n_out * dl) = 0 for a closed cell.
		error[elemL][0] += In[i][0]*length;
		error[elemL][1] += In[i][1]*length;
		error[elemR][0] -= In[i][0]*length;
		error[elemR][1] -= In[i][1]*length;
		
	}
	// loop through boundary edges
	for (int i = 0; i < B2E_size; i++){
		int Belem = B2E[i][0] - 1;
		int face = B2E[i][1];
		// int bgroup = B2E[i][2];
		double length = -1;
		vector<double> coord1(2);
		vector<double> coord2(2);

		// calculate length of face
		if (face == 1){
			int node1 = elem[Belem][1] - 1;
			int node2 = elem[Belem][2] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		else if (face == 2){
			int node1 = elem[Belem][0] - 1;
			int node2 = elem[Belem][2] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		else if (face == 3){
			int node1 = elem[Belem][0] - 1;
			int node2 = elem[Belem][1] - 1;
			coord1[0] = nodes[node1][0];
			coord2[0] = nodes[node2][0];
			coord1[1] = nodes[node1][1];
			coord2[1] = nodes[node2][1];

			length = sqrt(pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2));
		}
		
		// compute error of boundary edges with length and normal vectors
		error[Belem][0] += Bn[i][0]*length;
		error[Belem][1] += Bn[i][1]*length;

	}
	// compute magnitude of error for each element
	for (int i = 0; i < nelem; i++){
		err_mag[i] = sqrt(pow(error[i][0],2) + pow(error[i][1],2));
	}
	return err_mag;
}

#endif /* processMesh_h */