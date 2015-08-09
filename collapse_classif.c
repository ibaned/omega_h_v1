#include "collapse_classif.h"
#include "collapse_codes.h"
#include "tables.h"
#include "ints.h"
#include <assert.h>

void check_collapse_classif(
    unsigned elem_dim,
    unsigned nedges,
    unsigned* col_codes,
    unsigned const* class_dim_of_verts,
    unsigned const* verts_of_elems,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_verts_offsets,
    unsigned const* verts_of_verts,
    unsigned const* elems_of_edges_offsets,
    unsigned const* elems_of_edges,
    unsigned const* elems_of_edges_directions)
{
  assert(elem_dim >= 2);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned opp_dim = get_opposite_dim(elem_dim, 1);
  unsigned verts_per_opp = the_down_degrees[opp_dim][0];
  unsigned const* elem_ents_opp_edges =
    the_opposite_orders[elem_dim][1];
  unsigned const* const* elem_verts_of_opps =
    the_canonical_orders[elem_dim][opp_dim][0];
  for (unsigned i = 0; i < nedges; ++i) {
    if (col_codes[i] == DONT_COLLAPSE)
      continue;
    unsigned const* verts_of_edge = verts_of_edges + i * 2;
    unsigned class_dims_of_edge[2];
    class_dims_of_edge[0] = class_dim_of_verts[verts_of_edge[0]];
    class_dims_of_edge[1] = class_dim_of_verts[verts_of_edge[1]];
    if (class_dims_of_edge[0] == class_dims_of_edge[1])
      continue; /* equal orders, rest doesn't matter ? */
    for (unsigned j = 0; j < 2; ++j) {
      if (!(col_codes[i] & (1<<j)))
        continue;
      unsigned col_dim = class_dims_of_edge[j];
      unsigned gen_dim = class_dims_of_edge[1 - j];
      if (col_dim < gen_dim)
        col_codes[i] &= ~(1<<j);
    }
    if (col_codes[i] == DONT_COLLAPSE)
      continue;
    unsigned ring_buf[MAX_UP];
    unsigned ring_buf_size = 0;
    unsigned first_elem_use = elems_of_edges_offsets[i];
    unsigned end_elem_use = elems_of_edges_offsets[i + 1];
    for (unsigned j = first_elem_use; j < end_elem_use; ++j) {
      unsigned elem = elems_of_edges[j];
      unsigned direction = elems_of_edges_directions[j];
      unsigned opp = elem_ents_opp_edges[direction];
      unsigned const* elem_verts_of_opp = elem_verts_of_opps[opp];
      unsigned const* verts_of_elem = verts_of_elems + elem * verts_per_elem;
      for (unsigned k = 0; k < verts_per_opp; ++k) {
        unsigned vert = verts_of_elem[elem_verts_of_opp[k]];
        ring_buf_size = add_unique(ring_buf, ring_buf_size, vert);
      }
    }
    assert(ring_buf_size);
    for (unsigned j = 0; j < 2; ++j) {
      if (!(col_codes[i] & (1<<j)))
        continue;
      unsigned col_vert = verts_of_edge[j];
      unsigned gen_vert = verts_of_edge[1 - j];
      unsigned first_vert_use = verts_of_verts_offsets[col_vert];
      unsigned end_vert_use = verts_of_verts_offsets[col_vert + 1];
      unsigned col_class_dim = class_dims_of_edge[j];
      unsigned has_lesser = 0;
      unsigned has_equal = 0;
      for (unsigned k = first_vert_use; k < end_vert_use; ++k) {
        unsigned vert = verts_of_verts[k];
        if (vert == gen_vert || has(ring_buf, ring_buf_size, vert))
          continue;
        unsigned class_dim = class_dim_of_verts[vert];
        if (class_dim == col_class_dim)
          has_equal = 1;
        else if (class_dim < col_class_dim)
          has_lesser = 1;
      }
      unsigned is_ok = (has_equal && !has_lesser);
      if (!is_ok)
        col_codes[i] &= ~(1<<j);
    }
  }
}
