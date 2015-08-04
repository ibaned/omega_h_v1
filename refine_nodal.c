#include "refine_nodal.h"
#include "average.h"
#include "tables.h"
#include <stdlib.h>

double* refine_nodal(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned comps_per_vert,
    double const* field)
{
  double* out = malloc(sizeof(double) * comps_per_vert * nsrcs);
  unsigned verts_per_src = the_down_degrees[src_dim][0];
  for (unsigned i = 0; i < nsrcs; ++i) {
    if (gen_offset_of_srcs[i] == gen_offset_of_srcs[i + 1])
      continue;
    average_element_field(
        verts_per_src,
        verts_of_srcs + i * verts_per_src,
        field,
        comps_per_vert,
        out + gen_offset_of_srcs[i] * comps_per_vert);
  }
  return out;
}
