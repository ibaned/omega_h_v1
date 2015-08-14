#include "bridge_graph.h"
#include "ints.h"
#include "tables.h"
#include <assert.h>
#include <stdlib.h>

void bridge_graph_general(
    unsigned nverts,
    unsigned const adj_offsets[],
    unsigned const adj[],
    unsigned* nedges_out,
    unsigned** verts_of_edges_out,
    unsigned** directions_out)
{
  unsigned nhalf_edges = adj_offsets[nverts];
  assert(nhalf_edges % 2 == 0);
  unsigned nedges = nhalf_edges / 2;
  unsigned* degree_of_verts = malloc(sizeof(unsigned) * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_adj = adj_offsets[i];
    unsigned end_adj = adj_offsets[i + 1];
    unsigned degree_of_vert = 0;
    for (unsigned j = first_adj; j < end_adj; ++j)
      if (i < adj[j])
        ++degree_of_vert;
    degree_of_verts[i] = degree_of_vert;
  }
  unsigned* bridge_offsets = ints_exscan(degree_of_verts, nverts);
  free(degree_of_verts);
  unsigned* verts_of_edges = malloc(sizeof(unsigned) * nedges * 2);
  unsigned* directions = 0;
  if (directions_out)
    directions = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_adj = adj_offsets[i];
    unsigned end_adj = adj_offsets[i + 1];
    unsigned edge = bridge_offsets[i];
    for (unsigned j = first_adj; j < end_adj; ++j)
      if (i < adj[j]) {
        verts_of_edges[edge * 2 + 0] = i;
        verts_of_edges[edge * 2 + 1] = adj[j];
        ++edge;
        if (directions)
          directions[edge] = j - first_adj;
      }
  }
  free(bridge_offsets);
  *nedges_out = nedges;
  *verts_of_edges_out = verts_of_edges;
  if (directions_out)
    *directions_out = directions;
}

void bridge_graph(
    unsigned nverts,
    unsigned const adj_offsets[],
    unsigned const adj[],
    unsigned* nedges_out,
    unsigned** verts_of_edges_out)
{
  bridge_graph_general(nverts, adj_offsets, adj,
      nedges_out, verts_of_edges_out, 0);
}

void bridge_dual_graph(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const elems_to_elems[],
    unsigned* nfaces_out,
    unsigned** elems_of_faces_out,
    unsigned** elem_face_of_faces_out)
{
  unsigned faces_per_elem = the_down_degrees[elem_dim][elem_dim - 1];
  unsigned* degrees = malloc(sizeof(unsigned) * nelems);
  ints_fill(degrees, nelems, faces_per_elem);
  unsigned* offsets = ints_exscan(degrees, nelems);
  free(degrees);
  bridge_graph_general(nelems, offsets, elems_to_elems,
      nfaces_out, elem_face_of_faces_out, elem_face_of_faces_out);
  free(offsets);
}
