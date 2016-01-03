#include "swap_topology.h"

#include <assert.h>

#include "edge_ring.h"
#include "edge_swap.h"
#include "ints.h"
#include "loop.h"

unsigned* get_swap_topology_offsets(
    unsigned ent_dim,
    unsigned nedges,
    unsigned const* indset,
    unsigned const* ring_sizes)
{
  assert(ent_dim == 3);
  unsigned* nents_of_edge = LOOP_MALLOC(unsigned, nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    nents_of_edge[i] = indset[i] ?
      count_swap_ents(ring_sizes[i], ent_dim) : 0;
    if (indset[i]) {
      assert(nents_of_edge[i] >= 2);
      assert(nents_of_edge[i] <= 10);
    } else {
      assert(nents_of_edge[i] == 0);
    }
  }
  unsigned* offsets = uints_exscan(nents_of_edge, nedges);
  loop_free(nents_of_edge);
  return offsets;
}

unsigned* swap_topology(
    unsigned ent_dim,
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
  unsigned* out = LOOP_MALLOC(unsigned, ngen_elems * 4);
  for (unsigned i = 0; i < nedges; ++i) {
    if (!candidates[i])
      continue;
    unsigned* edge_out = out + gen_offset_of_edges[i] * 4;
    unsigned edge_v[2];
    unsigned ring_v[MAX_EDGE_SWAP+1];
    unsigned ring_size = find_edge_ring(i,
        tets_of_edges_offsets, tets_of_edges, tets_of_edges_directions,
        verts_of_edges, verts_of_tets, edge_v, ring_v);
    assert(ring_size <= MAX_EDGE_SWAP);
    unsigned ngen_elems_edge = 2 * swap_mesh_sizes[ring_size];
    assert(ngen_elems_edge == gen_offset_of_edges[i + 1] - gen_offset_of_edges[i]);
    get_swap_ents(ring_size, edge_codes[i], ent_dim, edge_v, ring_v, edge_out);
  }
  return out;
}
