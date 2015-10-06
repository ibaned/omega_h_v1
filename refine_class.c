#include "refine_class.h"
#include "tables.h"
#include "loop.h"

void refine_class(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned const* dims_in,
    unsigned const* ids_in,
    unsigned** p_dims_out,
    unsigned** p_ids_out)
{
  unsigned nsplits = gen_offset_of_srcs[nsrcs];
  unsigned* dims_out = loop_malloc(sizeof(unsigned) * nsplits);
  unsigned* ids_out = loop_malloc(sizeof(unsigned) * nsplits);
  unsigned verts_per_src = the_down_degrees[src_dim][0];
  for (unsigned i = 0; i < nsrcs; ++i) {
    if (gen_offset_of_srcs[i] == gen_offset_of_srcs[i + 1])
      continue;
    unsigned const* verts_of_src = verts_of_srcs + i * verts_per_src;
    unsigned dim = INVALID;
    unsigned id = INVALID;
    for (unsigned j = 0; j < verts_per_src; ++j) {
      unsigned vert = verts_of_src[j];
      unsigned vert_dim = dims_in[vert];
      if (j == 0 || vert_dim > dim) {
        dim = vert_dim;
        id = ids_in[vert];
      }
    }
    dims_out[gen_offset_of_srcs[i]] = dim;
    ids_out[gen_offset_of_srcs[i]] = id;
  }
  *p_dims_out = dims_out;
  *p_ids_out = ids_out;
}
