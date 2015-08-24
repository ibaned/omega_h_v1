#include "warp_to_limit.h"
#include "doubles.h"
#include "quality.h"
#include <stdlib.h>
#include <stdio.h>

void warp_to_limit(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords,
    double const* warps,
    double** p_coords,
    double** p_warps)
{
  double* coords_out = malloc(sizeof(double) * 3 * nverts);
  double factor = 1;
  double* warps_out = 0;
  doubles_axpy(factor, warps, coords, coords_out, 3 * nverts);
  while (min_element_quality(
        elem_dim, nelems, verts_of_elems, coords_out) <= 0) {
    factor /= 2;
    doubles_axpy(factor, warps, coords, coords_out, 3 * nverts);
  }
  *p_coords = coords_out;
  if (p_warps) {
    warps_out = malloc(sizeof(double) * 3 * nverts);
    doubles_axpy(-factor, warps, warps, warps_out, 3 * nverts);
    *p_warps = warps_out;
  }
}
