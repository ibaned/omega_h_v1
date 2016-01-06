#include "check_collapse_class.h"

#include <assert.h>

#include "collapse_codes.h"
#include "ints.h"
#include "mesh.h"
#include "tables.h"

static void check_reduced_collapse_class(struct mesh* m, unsigned* col_codes)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nedges = mesh_count(m, 1);
  unsigned const* class_dim_of_verts = mesh_find_tag(
      m, 0, "class_dim")->d.u32;
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* verts_of_verts_offsets =
      mesh_ask_star(m, 0, elem_dim)->offsets;
  unsigned const* verts_of_verts =
      mesh_ask_star(m, 0, elem_dim)->adj;
  unsigned const* elems_of_edges_offsets =
      mesh_ask_up(m, 1, elem_dim)->offsets;
  unsigned const* elems_of_edges =
      mesh_ask_up(m, 1, elem_dim)->adj;
  unsigned const* elems_of_edges_directions =
      mesh_ask_up(m, 1, elem_dim)->directions;
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
      if (!collapses(col_codes[i], j))
        continue;
      unsigned col_dim = class_dims_of_edge[j];
      unsigned gen_dim = class_dims_of_edge[1 - j];
      if (col_dim < gen_dim)
        col_codes[i] = dont_collapse(col_codes[i], j);
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
      if (!collapses(col_codes[i], j))
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
        col_codes[i] = dont_collapse(col_codes[i], j);
    }
  }
}

/* for now we'll just do the simple check and skip
   things like rings and exposed curved faces */
/* TODO: at least check the exposed face issue:

   |\ R1
   |R\
   |0/
   |/<-collapse */

static void check_full_collapse_class(struct mesh* m, unsigned* col_codes)
{
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* class_dim_of_verts = mesh_find_tag(
      m, 0, "class_dim")->d.u32;
  unsigned const* class_dim_of_edges = mesh_find_tag(
      m, 1, "class_dim")->d.u32;
  for (unsigned i = 0; i < nedges; ++i) {
    if (col_codes[i] == DONT_COLLAPSE)
      continue;
    unsigned const* verts_of_edge = verts_of_edges + i * 2;
    unsigned class_dims[3];
    enum { V0, V1, E };
    class_dims[V0] = class_dim_of_verts[verts_of_edge[0]];
    class_dims[V1] = class_dim_of_verts[verts_of_edge[1]];
    class_dims[E] = class_dim_of_edges[i];
    if (class_dims[V0] == class_dims[V1]) {
      if (class_dims[V0] != class_dims[E])
        col_codes[i] = DONT_COLLAPSE;
      continue;
    }
    for (unsigned j = 0; j < 2; ++j) {
      if (!collapses(col_codes[i], j))
        continue;
      if (class_dims[j] != class_dims[E])
        col_codes[i] = dont_collapse(col_codes[i], j);
    }
  }
}

void check_collapse_class(struct mesh* m, unsigned* col_codes)
{
  if (mesh_get_rep(m) == MESH_REDUCED)
    check_reduced_collapse_class(m, col_codes);
  else
    check_full_collapse_class(m, col_codes);
}
