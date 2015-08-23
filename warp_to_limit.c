#include "warp_to_limit.h"
#include "doubles.h"
#include "element_qualities.h"
#include <stdlib.h>
#include <stdio.h>

double* warp_to_limit(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords,
    double const* warps)
{
  double* out = malloc(sizeof(double) * 3 * nverts);
  double factor = 1;
  doubles_axpy(factor, warps, coords, out, 3 * nverts);
  while (min_element_quality(elem_dim, nelems, verts_of_elems, out) <= 0) {
    factor /= 2;
    doubles_axpy(factor, warps, coords, out, 3 * nverts);
  }
  printf("final factor %f\n", factor);
  return out;
}
