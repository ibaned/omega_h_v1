#include "bridge_graph.h"

#include <assert.h>

#include "ints.h"
#include "loop.h"
#include "tables.h"

static void bridge_graph_general(
    unsigned nverts,
    unsigned const* adj_offsets,
    unsigned const* adj,
    unsigned* nedges_out,
    unsigned** verts_of_edges_out,
    unsigned** directions_out)
{
  unsigned nhalf_edges = adj_offsets[nverts];
  assert(nhalf_edges % 2 == 0);
  unsigned* degree_of_verts = loop_malloc(sizeof(unsigned) * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_adj = adj_offsets[i];
    unsigned end_adj = adj_offsets[i + 1];
    unsigned degree_of_vert = 0;
    for (unsigned j = first_adj; j < end_adj; ++j)
      if (i < adj[j])
        ++degree_of_vert;
    degree_of_verts[i] = degree_of_vert;
  }
  unsigned* bridge_offsets = uints_exscan(degree_of_verts, nverts);
  loop_free(degree_of_verts);
  unsigned nedges = bridge_offsets[nverts];
  unsigned* verts_of_edges = loop_malloc(sizeof(unsigned) * nedges * 2);
  unsigned* directions = 0;
  if (directions_out)
    directions = loop_malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_adj = adj_offsets[i];
    unsigned end_adj = adj_offsets[i + 1];
    unsigned edge = bridge_offsets[i];
    for (unsigned j = first_adj; j < end_adj; ++j)
      if (i < adj[j]) {
        verts_of_edges[edge * 2 + 0] = i;
        verts_of_edges[edge * 2 + 1] = adj[j];
        if (directions)
          directions[edge] = j - first_adj;
        ++edge;
      }
  }
  loop_free(bridge_offsets);
  *nedges_out = nedges;
  *verts_of_edges_out = verts_of_edges;
  if (directions_out)
    *directions_out = directions;
}

void bridge_graph(
    unsigned nverts,
    unsigned const* adj_offsets,
    unsigned const* adj,
    unsigned* nedges_out,
    unsigned** verts_of_edges_out)
{
  bridge_graph_general(nverts, adj_offsets, adj,
      nedges_out, verts_of_edges_out, 0);
}

void bridge_dual_graph(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* elems_of_elems,
    unsigned* nfaces_out,
    unsigned** elems_of_faces_out,
    unsigned** elem_face_of_faces_out)
{
  unsigned faces_per_elem = the_down_degrees[elem_dim][elem_dim - 1];
  unsigned* degrees = loop_malloc(sizeof(unsigned) * nelems);
  uints_fill(degrees, nelems, faces_per_elem);
  unsigned* offsets = uints_exscan(degrees, nelems);
  loop_free(degrees);
  bridge_graph_general(nelems, offsets, elems_of_elems,
      nfaces_out, elems_of_faces_out, elem_face_of_faces_out);
  loop_free(offsets);
}
