#include "derive_edges.h"
#include "up_from_down.h"
#include "bridge_graph.h"
#include "star.h"
#include <stdlib.h>

void derive_edges(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    unsigned* nedges_out,
    unsigned** verts_of_edges_out)
{
  unsigned* elems_of_verts_offsets;
  unsigned* elems_of_verts;
  up_from_down(elem_dim, 0, nelems, nverts, verts_of_elems,
      &elems_of_verts_offsets, &elems_of_verts, 0);
  unsigned* verts_of_verts_offsets;
  unsigned* verts_of_verts;
  get_star(0, elem_dim, nverts, elems_of_verts_offsets, elems_of_verts,
      verts_of_elems, &verts_of_verts_offsets, &verts_of_verts);
  free(elems_of_verts_offsets);
  free(elems_of_verts);
  bridge_graph(nverts, verts_of_verts_offsets, verts_of_verts,
      nedges_out, verts_of_edges_out);
  free(verts_of_verts_offsets);
  free(verts_of_verts);
}
