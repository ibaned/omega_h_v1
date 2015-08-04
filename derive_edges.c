#include "derive_edges.h"
#include "up_from_down.h"
#include "star.h"
#include "bridge_graph.h"
#include <stdlib.h>

void derive_edges(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    unsigned* nedges_out,
    unsigned** verts_of_edges_out)
{
  struct up_adj vert_elems = up_from_down(
      elem_dim,
      0,
      nelems,
      nverts,
      verts_of_elems);
  free(vert_elems.directions); /* don't care about these */
  struct star vert_verts = get_star(
      0,
      elem_dim,
      nverts,
      vert_elems.offsets,
      vert_elems.edges,
      verts_of_elems);
  free(vert_elems.offsets);
  free(vert_elems.edges);
  bridge_graph(
      nverts,
      vert_verts.offsets,
      vert_verts.edges,
      nedges_out,
      verts_of_edges_out);
  free(vert_verts.offsets);
  free(vert_verts.edges);
}
