#include "swap_topology.h"
#include "edge_swap.h"
#include "edge_ring.h"
#include <stdlib.h>
#include <assert.h>

unsigned* swap_topology(
    unsigned nedges,
    unsigned const* candidates,
    unsigned const* gen_offset_of_edges,
    unsigned const* edge_codes,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets)
{
  unsigned ngen_elems = gen_offset_of_edges[nedges];
  unsigned* out = malloc(sizeof(unsigned) * ngen_elems * 4);
  for (unsigned i = 0; i < nedges; ++i) {
    if (!candidates[i])
      continue;
    unsigned* edge_out = out + gen_offset_of_edges[i];
    unsigned edge_v[2];
    unsigned ring_v[MAX_EDGE_SWAP];
    unsigned ring_size = find_edge_ring(i,
        tets_of_edges_offsets, tets_of_edges, tets_of_edges_directions,
        verts_of_edges, verts_of_tets, edge_v, ring_v);
    assert(ring_size <= MAX_EDGE_SWAP);
    unsigned ngen_elems_edge = 2 * swap_mesh_sizes[ring_size];
    assert(ngen_elems_edge = gen_offset_of_edges[i + 1] - gen_offset_of_edges[i]);
    apply_edge_swap(ring_size, edge_codes[i], edge_v, ring_v, edge_out);
  }
  return out;
}
