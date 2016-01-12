#include "bridge_graph.h"

#include <assert.h>

#include "ints.h"
#include "loop.h"
#include "tables.h"

LOOP_KERNEL( degree_count,
	unsigned const* adj,
    unsigned const* adj_offsets,
    unsigned* degree_of_verts)

  unsigned first_adj = adj_offsets[i];
  unsigned end_adj = adj_offsets[i+1];
  unsigned degree_of_vert = 0;
  for( unsigned j = first_adj ; j< end_adj ; j++) // Could make secondary Kernal Launch
    if( i < adj[j])
      ++degree_of_vert;
  degree_of_verts[i] = degree_of_vert;
}

LOOP_KERNEL( edge_count,
    unsigned const* adj,
    unsigned const* adj_offsets,
    unsigned* bridge_offsets,
    unsigned* verts_of_edges,
    unsigned* directions)

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


static void bridge_graph_general(
    unsigned nverts,
    unsigned const* adj_offsets,
    unsigned const* adj,
    unsigned* nedges_out,
    unsigned** verts_of_edges_out,
    unsigned** directions_out)
{
  unsigned* degree_of_verts = LOOP_MALLOC(unsigned, nverts);
  LOOP_EXEC(degree_count,
    nverts,
    adj,
    adj_offsets,
    degree_of_verts);
  unsigned* bridge_offsets = uints_exscan(degree_of_verts, nverts);
  loop_free(degree_of_verts);
  unsigned nedges = bridge_offsets[nverts];
  unsigned* verts_of_edges = LOOP_MALLOC(unsigned, nedges * 2);
  unsigned* directions = 0;
  if (directions_out)
    directions = LOOP_MALLOC(unsigned, nedges);
  LOOP_EXEC(edge_count,
    nverts,
    adj,
    adj_offsets,
    bridge_offsets,
    verts_of_edges,
    directions);
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
    unsigned* nsides_out,
    unsigned** elems_of_sides_out,
    unsigned** elem_side_of_sides_out)
{
  unsigned sides_per_elem = the_down_degrees[elem_dim][elem_dim - 1];
  unsigned* degrees = uints_filled(nelems, sides_per_elem);
  unsigned* offsets = uints_exscan(degrees, nelems);
  loop_free(degrees);
  bridge_graph_general(nelems, offsets, elems_of_elems,
      nsides_out, elems_of_sides_out, elem_side_of_sides_out);
  loop_free(offsets);
}
