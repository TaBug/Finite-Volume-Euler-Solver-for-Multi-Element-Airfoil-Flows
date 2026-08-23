#ifndef LIMITERS_H
#define LIMITERS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
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

#endif
