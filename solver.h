#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <Eigen/Dense>
#include "processMesh.h"
#include "fluxes.h"
using namespace std;
using namespace Eigen;

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

// Computes the length of an edge given the global edge index and the node coordinates
double computeEdgeLength(int iFace, vector<vector<double>> const& face, vector<vector<double>> const& nodes) {
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
	vector<double> alpha;  // [nelem][4], flattened - one limiter value per state variable
	explicit fvWorkspace(int nelem)
		: grad(size_t(nelem) * 8, 0.0), alpha(size_t(nelem) * 4, 1.0) {}
};

// State on the far side of a boundary face, written as four doubles.
//
// Wall: an inviscid wall allows slip but no through-flow, so the normal component
// of the velocity is projected out, v_b = v - (v.n)n. Density and total energy are
// carried over from the interior cell on purpose - holding rho*E fixed while the
// velocity shrinks is what lets the flux routine recover the wall pressure as
// p_b = (gamma-1)(rho*E - 0.5*rho*|v_b|^2). The normal must be a UNIT vector: the
// projection only removes the full normal component when |n| = 1, and silently
// under-corrects otherwise.
//
// Farfield: the full freestream state, valid while the outer boundary is far
// enough out that no characteristic information needs to travel back in.
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

// Included here rather than at the top of the file: the limiters need
// fvWorkspace complete, and the LCD limiter builds the same ghost state the
// residual does, so boundaryStateFlat has to be declared before it too.
#include "limiters.h"

// The flux is a template parameter, not a switch inside the loop: a function
// pointer known at compile time is a direct, inlinable call, so the flux body
// is folded into the face loop instead of costing a call per face.
template <double (*FluxCore)(const double*, const double*, double, double, double, double*)>
static void secondOrderResidualImpl(meshData const& mesh, meshGeom const& g,
									fvWorkspace& ws, const double* u, double* R, string const& limiterType) {
	const int nelem = mesh.nelem;
	double* grad = ws.grad.data();
	fill(ws.grad.begin(), ws.grad.end(), 0.0); // clear the gradient workspace
	fill(R, R + size_t(nelem) * 5, 0.0); // clear the residual workspace

	// Pass 1: grad(u)_i = (1/A_i) * sum over faces of uhat * n_out * dl
	//
	// Skipped for LCD: that limiter fits its own gradient - the plane through the
	// three neighbouring centroids - and overwrites every entry of ws.grad, so
	// computing Green-Gauss first would be thrown away in full.
	if (limiterType != "LCD") {
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
	}

	// One limiter value per state variable per element. Both limiters write every
	// entry themselves, so alpha only needs resetting when no limiter runs and the
	// reconstruction below must stay unlimited. computeLCDLimiter additionally
	// fills ws.grad, which is why the Green-Gauss pass above is skipped for it.
	double* alpha = ws.alpha.data();
	if (limiterType == "BJ") {
		computeBarthJespersenLimiter(alpha, mesh, u, ws);
	} 
	else if (limiterType == "LCD") {
		computeLCDLimiter(alpha, mesh, u, g, ws);
	}
	else {
		fill(ws.alpha.begin(), ws.alpha.end(), 1.0);
	}

	// Pass 2: reconstruct to each face midpoint and accumulate the flux residual
	for (size_t i = 0; i < g.faces.size(); i++) {
		faceCache const& f = g.faces[i];
		const double* gL = grad + size_t(f.iElemL) * 8;
		const double* cL = u + size_t(f.iElemL) * 4;
		const double* alphaL = alpha + size_t(f.iElemL) * 4;
		double uLr[4], uRr[4];
		for (int k = 0; k < 4; k++)
			uLr[k] = cL[k] + alphaL[k] * (gL[k * 2] * f.dxL + gL[k * 2 + 1] * f.dyL); // reconstructed left state at the face midpoint

		if (f.iElemR >= 0) {
			const double* gR = grad + size_t(f.iElemR) * 8;
			const double* cR = u + size_t(f.iElemR) * 4;
			const double* alphaR = alpha + size_t(f.iElemR) * 4;
			for (int k = 0; k < 4; k++)
				uRr[k] = cR[k] + alphaR[k] * (gR[k * 2] * f.dxR + gR[k * 2 + 1] * f.dyR);
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
	if (limiterType != "NONE" && limiterType != "BJ" && limiterType != "LCD") {
		cerr << "ERROR: limiterType \"" << limiterType << "\" is not implemented.\n"
			 << "       Available: \"NONE\" (unlimited), \"BJ\" (Barth-Jespersen), \"LCD\" (Limited Central Difference).\n";
		exit(EXIT_FAILURE);
	}
	switch (opt) {
		case FLUX_ROE:     secondOrderResidualImpl<roeCore>(mesh, g, ws, u, R, limiterType);     return;
		case FLUX_RUSANOV: secondOrderResidualImpl<rusanovCore>(mesh, g, ws, u, R, limiterType); return;
		case FLUX_HLLE:    secondOrderResidualImpl<hlleCore>(mesh, g, ws, u, R, limiterType);    return;
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

#endif
