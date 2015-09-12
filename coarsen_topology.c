#include "coarsen_topology.h"
#include "tables.h"
#include <stdlib.h>

void coarsen_topology(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* gen_offset_of_elems,
    unsigned const* gen_vert_of_elems,
    unsigned const* gen_direction_of_elems,
    unsigned* ngen_elems_out,
    unsigned** verts_of_gen_elems_out)
{
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned base_dim = get_opposite_dim(elem_dim, 0);
  unsigned ngen_elems = gen_offset_of_elems[nelems];
  unsigned const* const* elem_verts_of_bases =
    the_canonical_orders[elem_dim][base_dim][0];
  unsigned const* elem_bases_opp_verts =
    the_opposite_orders[elem_dim][0];
  unsigned* verts_of_gen_elems = malloc(
      sizeof(unsigned) * ngen_elems * verts_per_elem);
  for (unsigned i = 0; i < nelems; ++i) {
    if (gen_offset_of_elems[i] == gen_offset_of_elems[i + 1])
      continue;
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    unsigned* verts_of_gen_elem = verts_of_gen_elems +
      gen_offset_of_elems[i] * verts_per_elem;
    unsigned direction = gen_direction_of_elems[i];
    unsigned base = elem_bases_opp_verts[direction];
    unsigned const* elem_verts_of_base = elem_verts_of_bases[base];
    for (unsigned j = 0; j < (verts_per_elem - 1); ++j)
      verts_of_gen_elem[j] = verts_of_elem[elem_verts_of_base[j]];
    verts_of_gen_elem[verts_per_elem - 1] = gen_vert_of_elems[i];
  }
  *ngen_elems_out = ngen_elems;
  *verts_of_gen_elems_out = verts_of_gen_elems;
}
