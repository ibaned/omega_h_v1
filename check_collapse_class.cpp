#include "check_collapse_class.hpp"

#include <cassert>

#include "collapse_codes.hpp"
#include "ints.hpp"
#include "mesh.hpp"
#include "tables.hpp"

/* for now we'll just do the simple check and skip
   things like rings and exposed curved faces */
/* TODO: at least check the exposed face issue:

   |\ R1
   |R\
   |0/
   |/<-collapse */

LOOP_KERNEL(check_edge_collapse_class,
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

static void check_full_collapse_class(struct mesh* m, unsigned* col_codes)
{
  assert(mesh_get_rep(m) == MESH_FULL);
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* class_dim_of_verts = mesh_find_tag(
      m, 0, "class_dim")->d.u32;
  unsigned const* class_dim_of_edges = mesh_find_tag(
      m, 1, "class_dim")->d.u32;
  LOOP_EXEC(check_edge_collapse_class, nedges,
      verts_of_edges,
      class_dim_of_verts,
      class_dim_of_edges,
      col_codes);
}

void check_collapse_class(struct mesh* m, unsigned* col_codes)
{
  check_full_collapse_class(m, col_codes);
}
