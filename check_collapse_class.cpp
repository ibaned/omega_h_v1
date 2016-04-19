#include "check_collapse_class.hpp"

#include <cassert>

#include "collapse_codes.hpp"
#include "ints.hpp"
#include "mesh.hpp"
#include "tables.hpp"

namespace omega_h {

LOOP_KERNEL(check_basic_edge_class,
    unsigned const* verts_of_edges,
    unsigned const* class_dim_of_verts,
    unsigned const* class_dim_of_edges,
    unsigned* col_codes)

  if (col_codes[i] == DONT_COLLAPSE)
    return;
  unsigned const* verts_of_edge = verts_of_edges + i * 2;
  unsigned class_dims[3];
  enum { V0, V1, E };
  class_dims[V0] = class_dim_of_verts[verts_of_edge[0]];
  class_dims[V1] = class_dim_of_verts[verts_of_edge[1]];
  class_dims[E] = class_dim_of_edges[i];
  if (class_dims[V0] == class_dims[V1]) {
    if (class_dims[V0] != class_dims[E])
      col_codes[i] = DONT_COLLAPSE;
    return;
  }
  for (unsigned j = 0; j < 2; ++j) {
    if (!collapses(col_codes[i], j))
      continue;
    if (class_dims[j] != class_dims[E])
      col_codes[i] = dont_collapse(col_codes[i], j);
  }
}

static void check_basic_class(struct mesh* m, unsigned* col_codes)
{
  assert(mesh_get_rep(m) == MESH_FULL);
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* class_dim_of_verts = mesh_find_tag(
      m, 0, "class_dim")->d.u32;
  unsigned const* class_dim_of_edges = mesh_find_tag(
      m, 1, "class_dim")->d.u32;
  LOOP_EXEC(check_basic_edge_class, nedges,
      verts_of_edges,
      class_dim_of_verts,
      class_dim_of_edges,
      col_codes);
}

LOOP_KERNEL(check_for_exposed_side,
    unsigned ent_dim,
    unsigned const* ents_of_edges_offsets,
    unsigned const* ents_of_edges,
    unsigned const* ents_of_edges_directions,
    unsigned const* verts_of_ents,
    unsigned nverts_per_ent,
    unsigned const* verts_of_edges,
    unsigned const* sides_of_ents,
    unsigned nsides_per_ent,
    unsigned const* class_dim_of_ents,
    unsigned const* class_dim_of_sides,
    unsigned* col_codes)

  if (col_codes[i] == DONT_COLLAPSE)
    return;
  unsigned const* const* ent_verts_of_edges =
    the_canonical_orders[ent_dim][1][0];
  unsigned const* ent_sides_opp_verts =
    the_opposite_orders[ent_dim][0];
  unsigned f = ents_of_edges_offsets[i];
  unsigned e = ents_of_edges_offsets[i + 1];
  for (unsigned j = f; j < e; ++j) {
    unsigned ent = ents_of_edges[j];
    unsigned dir = ents_of_edges_directions[j];
    unsigned ent_class_dim = class_dim_of_ents[ent];
    unsigned const* ent_verts_of_edge =
      ent_verts_of_edges[dir];
    for (unsigned k = 0; k < 2; ++k) {
      unsigned eev = ent_verts_of_edge[k];
      unsigned va = verts_of_ents[ent * nverts_per_ent + eev];
      unsigned l;
      for (l = 0; l < 2; ++l) {
        unsigned vb = verts_of_edges[i * 2 + l];
        if (va == vb) {
          unsigned es = ent_sides_opp_verts[eev];
          unsigned side = sides_of_ents[ent * nsides_per_ent + es];
          unsigned side_class_dim = class_dim_of_sides[side];
          if (side_class_dim != ent_class_dim)
            col_codes[i] = dont_collapse(col_codes[i], l);
          break;
        }
      }
      assert(l != 2);
    }
  }
}

static void check_for_exposed_sides(struct mesh* m, unsigned ent_dim,
    unsigned* col_codes)
{
  unsigned const* ents_of_edges_offsets =
    mesh_ask_up(m, 1, ent_dim)->offsets;
  unsigned const* ents_of_edges =
    mesh_ask_up(m, 1, ent_dim)->adj;
  unsigned const* ents_of_edges_directions =
    mesh_ask_up(m, 1, ent_dim)->directions;
  unsigned const* verts_of_ents =
    mesh_ask_down(m, ent_dim, 0);
  unsigned nverts_per_ent = the_down_degrees[ent_dim][0];
  unsigned const* verts_of_edges =
    mesh_ask_down(m, 1, 0);
  unsigned const* sides_of_ents =
    mesh_ask_down(m, ent_dim, ent_dim - 1);
  unsigned nsides_per_ent = the_down_degrees[ent_dim][ent_dim - 1];
  unsigned const* class_dim_of_ents = mesh_find_tag(
      m, ent_dim, "class_dim")->d.u32;
  unsigned const* class_dim_of_sides = mesh_find_tag(
      m, ent_dim - 1, "class_dim")->d.u32;
  unsigned nedges = mesh_count(m, 1);
  LOOP_EXEC(check_for_exposed_side, nedges,
      ent_dim,
      ents_of_edges_offsets,
      ents_of_edges,
      ents_of_edges_directions,
      verts_of_ents,
      nverts_per_ent,
      verts_of_edges,
      sides_of_ents,
      nsides_per_ent,
      class_dim_of_ents,
      class_dim_of_sides,
      col_codes);
}

void check_collapse_class(struct mesh* m, unsigned* col_codes)
{
  check_basic_class(m, col_codes);
  for (unsigned d = 2; d <= mesh_dim(m); ++d)
    check_for_exposed_sides(m, d, col_codes);
}

}
