#include "refine_nodal.h"
#include "average.h"
#include "tables.h"
#include <stdlib.h>

double* refine_nodal(
    unsigned split_dim,
    unsigned nsplit,
    unsigned const* split_verts,
    unsigned const* split_offsets,
    unsigned ncomp,
    double const* field)
{
  double* out = malloc(sizeof(double) * ncomp * nsplit);
  unsigned nsplit_vert = the_down_degrees[split_dim][0];
  for (unsigned i = 0; i < nsplit; ++i) {
    if (split_offsets[i] == split_offsets[i + 1])
      continue;
    average_element_field(
        nsplit_vert,
        split_verts,
        field,
        ncomp,
        out + split_offsets[i] * ncomp);
  }
  return out;
}
