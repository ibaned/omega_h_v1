#include "refine_topology.h"
#include "tables.h"
#include <stdlib.h>
#include <assert.h>

struct refined_topology refine_topology(
    unsigned elem_dim,
    unsigned split_dim,
    unsigned ent_dim,
    unsigned nelem,
    unsigned const* elem_verts,
    unsigned const* elem_split_offset,
    unsigned const* elem_split_vert,
    unsigned const* elem_split_direction)
{
  struct refined_topology out = {0,0,0};
  assert(elem_dim <= 3);
  assert(elem_dim >= split_dim);
  assert(elem_dim >= ent_dim);
  assert(ent_dim > 0);
  unsigned base_dim = ent_dim - 1;
  unsigned opposite_dim = get_opposite_dim(elem_dim, base_dim);
  if (split_dim < opposite_dim)
    return out;
  unsigned split_degree = the_down_degrees[split_dim][opposite_dim];
  assert(split_degree);
  unsigned nselem = elem_split_offset[nelem];
  if (!nselem)
    return out;
  out.nent = nselem * split_degree;
  unsigned nent_vert = the_down_degrees[ent_dim][0];
  unsigned nelem_vert = the_down_degrees[elem_dim][0];
  out.ent_verts = malloc(sizeof(unsigned) * nselem * nent_vert);
  unsigned const* const* any_split_opposite =
    the_canonical_orders[elem_dim][split_dim][opposite_dim];
  unsigned const* const* any_base_vert =
    the_canonical_orders[elem_dim][base_dim][0];
  unsigned const* opposite_base = the_opposite_orders[elem_dim][opposite_dim];
  for (unsigned i = 0; i < nelem; ++i) {
    if (elem_split_offset[i] == elem_split_offset[i + 1])
      continue;
    unsigned direction = elem_split_direction[i];
    unsigned const* split_opposite = any_split_opposite[direction];
    unsigned const* elem_vert = elem_verts + i * nelem_vert;
    unsigned split_vert = elem_split_vert[i];
    unsigned* ent_vert = out.ent_verts + elem_split_offset[i] * nent_vert;
    for (unsigned j = 0; j < split_degree; ++j) {
      unsigned opposite = split_opposite[j];
      unsigned base = opposite_base[opposite];
      unsigned const* base_vert = any_base_vert[base];
      for (unsigned k = 0; k < nelem_vert; ++k)
        ent_vert[k] = elem_vert[base_vert[k]];
      ent_vert[nelem_vert] = split_vert;
      ent_vert += nent_vert;
    }
  }
  return out;
}
