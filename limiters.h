#ifndef LIMITERS_H
#define LIMITERS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
// #include <Eigen/Geometry>
#include "processMesh.h"
using namespace std;

// Barth-Jespersen limiter, vertex based.
//
// For each element: collect the min/max cell average over the element itself and
// its face neighbours, reconstruct to each of the three vertices with the
// already-area-divided gradient in ws.grad, and shrink alpha until no vertex value
// overshoots that range. One alpha per state variable per element.
//
// alpha must have room for nelem*4 doubles; it is fully written here (seeded to 1
// per element), so the caller does not have to pre-fill it.
//
// mesh.r_node_centroid holds node-minus-centroid offsets - see genCentroid2NodesVec.
inline void computeBarthJespersenLimiter(double* alpha, meshData const& mesh, const double* u, fvWorkspace& ws) {
	const double* grad = ws.grad.data();

	for (int iElem = 0; iElem < mesh.nelem; iElem++) {
		const double* u0 = u + size_t(iElem) * 4; // cell average of the current element
		const double* grad0 = grad + size_t(iElem) * 8; // gradient of the current element

		double* alphaElem = alpha + size_t(iElem) * 4; // limiter values for the current element

		// Initialize uMin and uMax with the cell average of the current element, and
		// alpha with 1 (unlimited) so the mins below can only tighten it.
		double uMin[4], uMax[4];
		for (int k = 0; k < 4; k++) {
			uMin[k] = u0[k];
			uMax[k] = u0[k];
			alphaElem[k] = 1.0;
		}

		for (int iFaceLocal = 0; iFaceLocal < 3; iFaceLocal++) {
			int iFaceGlobal = mesh.E2F[iElem][iFaceLocal]; // Get the global face index for the current local face

			// Global face indices run interior first, boundary second (see genE2F), and
			// I2E only has niedge rows - indexing it with a boundary face reads off the
			// end. A boundary face has no neighbouring cell, so it contributes nothing to
			// the min/max stencil; the element's own average keeps alpha well defined.
			if (iFaceGlobal >= mesh.niedge) continue;

			// I2E stores 1-based element numbers, iElem is 0-based - convert before comparing
			int elemL = int(mesh.I2E[iFaceGlobal][0]) - 1;
			int elemR = int(mesh.I2E[iFaceGlobal][2]) - 1;
			int elem_i = (elemL == iElem) ? elemR : elemL; // Get the neighboring element index
			const double* u_i = u + size_t(elem_i) * 4; // cell average of the neighboring element

			// Update uMin and uMax based on the neighboring element's cell average
			for (int k = 0; k < 4; k++) {
				uMin[k] = min(uMin[k], u_i[k]);
				uMax[k] = max(uMax[k], u_i[k]);
			}
		}

		for (int iNode = 0; iNode < 3; iNode++) {
			double rx = mesh.r_node_centroid[iElem][2 * iNode];
			double ry = mesh.r_node_centroid[iElem][2 * iNode + 1];

			for (int k = 0; k < 4; k++) {
				double uN_i_k = u0[k] + (rx * grad0[2 * k] + ry * grad0[2 * k + 1]); // Reconstructed value at the node for each state variable
				double d = uN_i_k - u0[k];
				
				// Roundoff floor. Over a patch where every average is bit-identical,
				// uMax - u0 is exactly 0 while the Green-Gauss sum still leaves a ~1e-16
				// gradient, so the unguarded ratio sets alpha to 0 on a perfectly smooth
				// region. The threshold is at the noise level, so a real overshoot -
				// which is orders of magnitude larger - is still limited exactly as before.
				double tol = 1e-12 * fabs(u0[k]) + 1e-300;
				
				if (d > tol) {
					alphaElem[k] = min(alphaElem[k], (uMax[k] - u0[k]) / d);
				}
				else if (d < -tol) {
					alphaElem[k] = min(alphaElem[k], (uMin[k] - u0[k]) / d);
				}
			}
		}

	}
}

// Limited Central Difference limiter, in the projected form of Hubbard (1999).
//
// Step (a), Hubbard Eqn 2.9: fit a plane through the three (centroid_x,
// centroid_y, u) points of the element's face neighbours - one plane per state
// variable - and take its slope as the candidate gradient L. This reproduces a
// linear field exactly.
//
// Step (b), Hubbard Eqn 2.12: the maximum-principle region is the set of
// gradients that create no new extremum at any face midpoint,
//
//     min(u_k - u_0, 0) <= r_0k . L <= max(u_k - u_0, 0),   k = 1,2,3
//
// where r_0k runs from this element's centroid to the midpoint of the edge it
// shares with neighbour k. In gradient space that is the intersection of three
// slabs: a convex polygon that always contains the origin, since the lower
// bound is never positive and the upper bound never negative.
//
// The limited gradient is THE POINT OF THAT REGION CLOSEST TO L - the
// projection variant Hubbard describes alongside Eqn 2.15. The plain LCD scheme
// instead rescales L along its own direction, which is cheaper but, as Hubbard
// notes in 2.1.1, "does not depend continuously on the solution data ... it may
// interfere with convergence to a steady state by causing limit cycling". The
// projection is continuous in the data and he reports it is more accurate too.
//
// Because the projected gradient is not parallel to L, no scalar can express
// it: this writes the limited gradient straight into ws.grad and leaves alpha
// at 1. The Green-Gauss gradient the residual would otherwise build is not the
// operator being limited here, so the caller skips that pass entirely and every
// element must be written, the degenerate case below included.
//
// A boundary face has no neighbouring cell, so it contributes the mirrored
// ghost point (centroid reflected through the face midpoint, carrying the same
// ghost state the residual builds) exactly as the reference computeP_boundary
// does. Hubbard instead drops boundary faces from both the fit and the extremum
// search; the ghost keeps this in step with the legacy solver.

// Euclidean projection of the gradient (gx, gy) onto the maximum-principle
// region { L : lo[k] <= r[k] . L <= hi[k] }, in place.
//
// The region is a convex polygon, so the closest point is either the gradient
// itself, the foot of the perpendicular onto one bounding line, or a corner
// where two of them meet. All three kinds are enumerated and the nearest
// feasible one wins; the origin seeds the search so a feasible answer always
// exists even when the region collapses to a point (a local extremum, where
// every u_k - u_0 shares a sign and the scheme drops to first order).
//
// tol is in the units of u: it absorbs the roundoff in r.L over a patch of
// bit-identical averages, where the region is exactly a point and an unguarded
// test would zero a perfectly smooth gradient.
inline void projectOntoMPRegion(const double rx[3], const double ry[3],
								const double lo[3], const double hi[3],
								double tol, double& gx, double& gy) {
	// Six half-planes a.L <= b: the three upper bounds, then the three lower
	// bounds negated into the same form.
	double ax[6], ay[6], b[6];
	for (int k = 0; k < 3; k++) {
		ax[k] = rx[k];      ay[k] = ry[k];      b[k] = hi[k];
		ax[k + 3] = -rx[k]; ay[k + 3] = -ry[k]; b[k + 3] = -lo[k];
	}

	// A NaN coordinate would pass every > comparison below and be mistaken for a
	// feasible point, so finiteness is part of feasibility.
	auto feasible = [&](double x, double y) {
		if (!isfinite(x) || !isfinite(y)) return false;
		for (int i = 0; i < 6; i++) {
			if (ax[i] * x + ay[i] * y > b[i] + tol) return false;
		}
		return true;
	};

	if (feasible(gx, gy)) return; // already inside - the common case, and no work

	double bestX = 0.0, bestY = 0.0;           // the origin is always feasible
	double bestD = gx * gx + gy * gy;          // ... at this distance

	// Candidates on the interior of an edge.
	for (int i = 0; i < 6; i++) {
		double n2 = ax[i] * ax[i] + ay[i] * ay[i];
		if (!(n2 > 0.0)) continue; // a zero-length offset constrains nothing
		double t = (ax[i] * gx + ay[i] * gy - b[i]) / n2;
		double px = gx - t * ax[i], py = gy - t * ay[i];
		if (!feasible(px, py)) continue;
		double d = (px - gx) * (px - gx) + (py - gy) * (py - gy);
		if (d < bestD) { bestD = d; bestX = px; bestY = py; }
	}

	// Candidates at a corner, where two bounding lines cross.
	for (int i = 0; i < 6; i++) {
		for (int j = i + 1; j < 6; j++) {
			double det = ax[i] * ay[j] - ay[i] * ax[j];
			// Near-parallel lines - which includes the two faces of one slab -
			// have no usable intersection. The test is on the sine of the angle
			// between them so it does not depend on how the rows are scaled.
			double scale = sqrt((ax[i] * ax[i] + ay[i] * ay[i]) * (ax[j] * ax[j] + ay[j] * ay[j]));
			if (!(fabs(det) > 1e-14 * scale)) continue;
			double px = (b[i] * ay[j] - ay[i] * b[j]) / det;
			double py = (ax[i] * b[j] - b[i] * ax[j]) / det;
			if (!feasible(px, py)) continue;
			double d = (px - gx) * (px - gx) + (py - gy) * (py - gy);
			if (d < bestD) { bestD = d; bestX = px; bestY = py; }
		}
	}

	gx = bestX;
	gy = bestY;
}
 
inline void computeLCDLimiter(double* alpha, meshData const& mesh, const double* u, meshGeom const& g,
							  fvWorkspace& ws) {
	double* grad = ws.grad.data();

	for (int iElem = 0; iElem < mesh.nelem; iElem++) {
		const double* u0 = u + size_t(iElem) * 4; // cell average of the current element

		double Ps[3][6];       // per local face: fit point x, y, then the four state values
		double rX[3], rY[3];   // per local face: centroid -> face midpoint, Hubbard's r_0k
		double uMaxFace[3][4]; // per local face and state: largest admissible increment (>= 0)
		double uMinFace[3][4]; // per local face and state: smallest admissible increment (<= 0)

		for (int iFaceLocal = 0; iFaceLocal < 3; iFaceLocal++) {
			int iFaceGlobal = mesh.E2F[iElem][iFaceLocal]; // Get the global face index for the current local face
			faceCache const& f = g.faces[iFaceGlobal];

			// dxL/dyL is the face midpoint measured from the LEFT element's
			// centroid; this element is the right one on roughly half its faces,
			// so pick the offset that is actually rooted at iElem.
			bool isLeft = (f.iElemL == iElem);
			double dx = isLeft ? f.dxL : f.dxR;
			double dy = isLeft ? f.dyL : f.dyR;
			rX[iFaceLocal] = dx;
			rY[iFaceLocal] = dy;

			// Fit point contributed by this face: the neighbour's centroid and
			// average, or the mirrored ghost when there is no neighbour.
			double ub[4];
			const double* u_i;
			if (f.iElemR < 0) { // boundary face - iElem is necessarily the left element
				boundaryStateFlat(u0, f.isWall != 0, g.uinf, f.nx, f.ny, ub);
				u_i = ub;
				// Ghost centroid = midpoint + (midpoint - centroid) = centroid + 2*(dx, dy)
				Ps[iFaceLocal][0] = mesh.centroids[iElem][0] + 2.0 * dx;
				Ps[iFaceLocal][1] = mesh.centroids[iElem][1] + 2.0 * dy;
			}
			else {
				int elem_i = isLeft ? f.iElemR : f.iElemL; // Get the neighboring element index
				u_i = u + size_t(elem_i) * 4; // cell average of the neighboring element
				Ps[iFaceLocal][0] = mesh.centroids[elem_i][0];
				Ps[iFaceLocal][1] = mesh.centroids[elem_i][1];
			}

			for (int k = 0; k < 4; k++) {
				Ps[iFaceLocal][k + 2] = u_i[k];
				uMaxFace[iFaceLocal][k] = max(u_i[k] - u0[k], 0.0);
				uMinFace[iFaceLocal][k] = min(u_i[k] - u0[k], 0.0);
			}
		}

		double* alphaElem = alpha + size_t(iElem) * 4;
		double* gradElem = grad + size_t(iElem) * 8;

		// The projected gradient carries all of the limiting, so alpha stays 1
		// and the residual's alpha * grad . r reduces to the limited gradient
		// itself. ws.alpha persists across RK stages, so this must be written
		// rather than assumed.
		for (int k = 0; k < 4; k++) alphaElem[k] = 1.0;

		// The in-plane part of the two fit edges does not depend on the state
		// variable, and neither does the z component of their cross product -
		// which is twice the signed area of the triangle the three fit points
		// span in (x, y), i.e. the common denominator of all four gradients.
		double v1x = Ps[1][0] - Ps[0][0], v1y = Ps[1][1] - Ps[0][1];
		double v2x = Ps[2][0] - Ps[0][0], v2y = Ps[2][1] - Ps[0][1];
		double nz = (v1x * v2y - v1y * v2x);
		
		// if (nz > 100) {
		// 	cout << "nz: " << nz << "\n";
		// }

		for (int k = 0; k < 4; k++) {
			double v1z = Ps[1][k + 2] - Ps[0][k + 2];
			double v2z = Ps[2][k + 2] - Ps[0][k + 2];
			// n = v1 x v2; the plane's gradient is (-nx/nz, -ny/nz)
			double nx = (v1y * v2z - v1z * v2y);
			double ny = (v1z * v2x - v1x * v2z);
			
			// The plane's gradient is (-nx/nz, -ny/nz), 
			// but if nz is very small the three points are nearly collinear and the plane is ill-defined.
			// In that case, we set the gradient to zero to avoid numerical issues.
			bool degenerate = fabs(nz) < 1e-10;
			double gx = degenerate ? 0.0 : -nx / nz;
			double gy = degenerate ? 0.0 : -ny / nz;

			double lo[3] = { uMinFace[0][k], uMinFace[1][k], uMinFace[2][k] };
			double hi[3] = { uMaxFace[0][k], uMaxFace[1][k], uMaxFace[2][k] };
			double tol = 1e-12 * fabs(u0[k]) + 1e-300;

			projectOntoMPRegion(rX, rY, lo, hi, tol, gx, gy);

			gradElem[2 * k] = gx;
			gradElem[2 * k + 1] = gy;
		}

		// for (int iFaceLocal = 0; iFaceLocal < 3; iFaceLocal++) {
		// 	int iFaceGlobal = mesh.E2F[iElem][iFaceLocal]; // Get the global face index for the current local face
		// 	faceCache const& f = g.faces[iFaceGlobal];

		// 	// dxL/dyL is the face midpoint measured from the LEFT element's
		// 	// centroid; this element is the right one on roughly half its faces,
		// 	// so pick the offset that is actually rooted at iElem.
		// 	bool isLeft = (f.iElemL == iElem);
		// 	double dx = isLeft ? f.dxL : f.dxR;
		// 	double dy = isLeft ? f.dyL : f.dyR;

		// 	for (int k = 0; k < 4; k++) {
		// 		double rL_k = dx * gradElem[2 * k] + dy * gradElem[2 * k + 1];

		// 		if (rL_k > uMaxFace[iFaceLocal][k]) {
		// 			alphaElem[k] = min(alphaElem[k], (uMaxFace[iFaceLocal][k]) / rL_k);
		// 		}
		// 		else if (rL_k < uMinFace[iFaceLocal][k]) {
		// 			alphaElem[k] = min(alphaElem[k], (uMinFace[iFaceLocal][k]) / rL_k);
		// 		}
		// 		else {
		// 			alphaElem[k] = 1.0; // No limiting needed for this face and state variable
		// 		}
		// 	}
		// }

	}
}
#endif
