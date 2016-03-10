#include "refine_nodal.hpp"

#include "arrays.hpp"
#include "element_field.hpp"
#include "loop.hpp"
#include "tables.hpp"

LOOP_KERNEL(refine_node,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned comps_per_vert,
    double const* field,
    unsigned verts_per_src,
    double* out)
  if (gen_offset_of_srcs[i] == gen_offset_of_srcs[i + 1])
    return;
  average_element_field(
      verts_per_src,
      verts_of_srcs + i * verts_per_src,
      field,
      comps_per_vert,
      out + gen_offset_of_srcs[i] * comps_per_vert);
}

double* refine_nodal(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned comps_per_vert,
    double const* field)
{
  unsigned nsplits = uints_at(gen_offset_of_srcs, nsrcs);
  double* out = LOOP_MALLOC(double, comps_per_vert * nsplits);
  unsigned verts_per_src = the_down_degrees[src_dim][0];
  LOOP_EXEC(refine_node, nsrcs,
      verts_of_srcs,
      gen_offset_of_srcs,
      comps_per_vert,
      field,
      verts_per_src,
      out);
  return out;
}
