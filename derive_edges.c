#include "derive_edges.h"
#include "up_from_down.h"
#include "star.h"
#include "bridge.h"
#include <stdlib.h>

struct derived_edges derive_edges(
    unsigned elem_dim,
    unsigned nelem,
    unsigned nvert,
    unsigned const* elem_verts)
{
  struct up_adj vert_elems = up_from_down(
      elem_dim,
      0,
      nelem,
      nvert,
      elem_verts);
  free(vert_elems.directions); /* don't care about these */
  struct star vert_verts = get_star(
      0,
      elem_dim,
      nvert,
      vert_elems.offsets,
      vert_elems.edges,
      elem_verts);
  free(vert_elems.offsets);
  free(vert_elems.edges);
  struct bridged_graph bg = bridge_graph(
      nvert,
      vert_verts.offsets,
      vert_verts.edges);
  free(vert_verts.offsets);
  free(vert_verts.edges);
  return (struct derived_edges) {
    .nedge = bg.nedge,
    .edge_verts = bg.edge_verts
  };
}
