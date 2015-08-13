#include "refine_classif.h"
#include "tables.h"
#include <stdlib.h>

unsigned* refine_classif(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned const* class_dim_of_verts)
{
  unsigned nsplits = gen_offset_of_srcs[nsrcs];
  unsigned* out = malloc(sizeof(double) * nsplits);
  unsigned verts_per_src = the_down_degrees[src_dim][0];
  for (unsigned i = 0; i < nsrcs; ++i) {
    if (gen_offset_of_srcs[i] == gen_offset_of_srcs[i + 1])
      continue;
    unsigned const* verts_of_src = verts_of_srcs + i * verts_per_src;
    unsigned class_dim = 4;
    for (unsigned j = 0; j < verts_per_src; ++j) {
      unsigned vert = verts_of_src[j];
      unsigned class_dim_of_vert = class_dim_of_verts[vert];
      if (class_dim_of_vert < class_dim)
        class_dim = class_dim_of_vert;
    }
    out[gen_offset_of_srcs[i]] = class_dim;
  }
  return out;
}
