#include "coarsen_topology.hpp"

#include "arrays.hpp"
#include "loop.hpp"
#include "tables.hpp"

LOOP_KERNEL(coarsen_elem,
    unsigned elem_dim,
    unsigned base_dim,
    unsigned verts_per_elem,
    unsigned const* gen_offset_of_elems,
    unsigned const* gen_vert_of_elems,
    unsigned const* gen_direction_of_elems,
    unsigned const* verts_of_elems,
    unsigned* verts_of_gen_elems)

  if (gen_offset_of_elems[i] == gen_offset_of_elems[i + 1])
    return;
  unsigned const* const* elem_verts_of_bases =
      the_canonical_orders[elem_dim][base_dim][0];
  unsigned const* elem_bases_opp_verts =
    the_opposite_orders[elem_dim][0];
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
  unsigned ngen_elems = array_at(gen_offset_of_elems, nelems);
  unsigned* verts_of_gen_elems = LOOP_MALLOC(unsigned,
      ngen_elems * verts_per_elem);
  LOOP_EXEC(coarsen_elem, nelems,
      elem_dim,
      base_dim,
      verts_per_elem,
      gen_offset_of_elems,
      gen_vert_of_elems,
      gen_direction_of_elems,
      verts_of_elems,
      verts_of_gen_elems);
  *ngen_elems_out = ngen_elems;
  *verts_of_gen_elems_out = verts_of_gen_elems;
}
