#include "edge_ring.h"
#include "edge_swap.h"
#include "tables.h"
#include <assert.h>

struct ev { unsigned a; unsigned b; };

unsigned find_edge_ring(
    unsigned edge,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    unsigned edge_v[2],
    unsigned ring_v[])
{
  unsigned first_use = tets_of_edges_offsets[edge];
  unsigned end_use = tets_of_edges_offsets[edge + 1];
  unsigned ring_size = end_use - first_use;
  assert(ring_size >= 3);
  if (ring_size > MAX_EDGE_SWAP)
    return ring_size;
  unsigned const* tet_edge_opp_edges = the_opposite_orders[3][1];
  unsigned const* const* tet_verts_of_edges = the_canonical_orders[3][1][0];
  edge_v[0] = verts_of_edges[edge * 2 + 0];
  edge_v[1] = verts_of_edges[edge * 2 + 1];
  struct ev tmp_ring[MAX_EDGE_SWAP];
  for (unsigned i = first_use; i < end_use; ++i) {
    unsigned tet = tets_of_edges[i];
    unsigned const* verts_of_tet = verts_of_tets + tet * 4;
    unsigned tet_edge = tets_of_edges_directions[i];
    unsigned const* tet_verts_of_edge = tet_verts_of_edges[tet_edge];
    unsigned align = 0;
    if (verts_of_tet[tet_verts_of_edge[0]] != edge_v[0])
      align = 1;
    assert(edge_v[align] = verts_of_tet[tet_verts_of_edge[0]]);
    unsigned tet_edge_opp = tet_edge_opp_edges[tet_edge];
    unsigned const* tet_verts_of_opp = tet_verts_of_edges[tet_edge_opp];
    tmp_ring[i - first_use].a = verts_of_tet[tet_verts_of_opp[align ^ 0]];
    tmp_ring[i - first_use].b = verts_of_tet[tet_verts_of_opp[align ^ 1]];
  }
  /* upward adjacency from edges to tets is not ordered
   * clockwise around, it is essentially random.
   * the following code uses insertion sort to
   * order the edges around the ring by matching their
   * endpoints.
   * this might actually be faster than doing something
   * with face adjacencies.
   * remember, at most 7 entries, no fancy sort needed
   */
  for (unsigned i = 0; i < ring_size - 1; ++i) {
    unsigned j;
    for (j = i + 1; j < ring_size; ++j)
      if (tmp_ring[j].a == tmp_ring[i].b)
        break;
    assert(j < ring_size);
    struct ev tmp = tmp_ring[j];
    tmp_ring[j] = tmp_ring[i + 1];
    tmp_ring[i + 1] = tmp;
  }
  for (unsigned i = 0; i < ring_size; ++i)
    ring_v[i] = tmp_ring[i].a;
  return ring_size;
}
