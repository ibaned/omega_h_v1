#include "swap_topology.hpp"

#include <cassert>

#include "algebra.hpp"
#include "arrays.hpp"
#include "edge_ring.hpp"
#include "edge_swap.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tables.hpp"

LOOP_KERNEL(edge_topo_offset,
    unsigned ent_dim,
    unsigned const* indset,
    unsigned const* ring_sizes,
    unsigned* nents_of_edges)
  nents_of_edges[i] = indset[i] ?
    count_swap_ents(ring_sizes[i], ent_dim) : 0;
}

unsigned* get_swap_topology_offsets(
    unsigned ent_dim,
    unsigned nedges,
    unsigned const* indset,
    unsigned const* ring_sizes)
{
  unsigned* nents_of_edges = LOOP_MALLOC(unsigned, nedges);
  LOOP_EXEC(edge_topo_offset, nedges,
      ent_dim,
      indset,
      ring_sizes,
      nents_of_edges);
  unsigned* offsets = uints_exscan(nents_of_edges, nedges);
  loop_free(nents_of_edges);
  return offsets;
}

LOOP_KERNEL(swap_cavity,
    unsigned ent_dim,
    unsigned const* candidates,
    unsigned const* gen_offset_of_edges,
    unsigned verts_per_ent,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    double const* coords,
    unsigned* out)
  if (!candidates[i])
    return;
  unsigned* edge_out = out + gen_offset_of_edges[i] * verts_per_ent;
  unsigned edge_v[2];
  unsigned ring_v[MAX_EDGE_SWAP+1];
  unsigned ring_size = find_edge_ring(i,
      tets_of_edges_offsets, tets_of_edges, tets_of_edges_directions,
      verts_of_edges, verts_of_tets, edge_v, ring_v);
  assert(ring_size <= MAX_EDGE_SWAP);
  double edge_x[2][3];
  copy_vector(coords + edge_v[0] * 3, edge_x[0], 3);
  copy_vector(coords + edge_v[1] * 3, edge_x[1], 3);
  double ring_x[MAX_EDGE_SWAP][3];
  for (unsigned j = 0; j < ring_size; ++j)
    copy_vector(coords + ring_v[j] * 3, ring_x[j], 3);
  struct swap_choice sc = choose_edge_swap(ring_size, edge_x, ring_x);
  get_swap_ents(ring_size, sc.code, ent_dim, edge_v, ring_v, edge_out);
}

static unsigned* swap_topology(
    unsigned ent_dim,
    unsigned nedges,
    unsigned const* candidates,
    unsigned const* gen_offset_of_edges,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    double const* coords)
{
  unsigned ngen_ents = uints_at(gen_offset_of_edges, nedges);
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned* out = LOOP_MALLOC(unsigned, ngen_ents * verts_per_ent);
  LOOP_EXEC(swap_cavity, nedges,
      ent_dim,
      candidates,
      gen_offset_of_edges,
      verts_per_ent,
      tets_of_edges_offsets,
      tets_of_edges,
      tets_of_edges_directions,
      verts_of_edges,
      verts_of_tets,
      coords,
      out);
  return out;
}

unsigned* mesh_swap_topology(
    struct mesh* m,
    unsigned ent_dim,
    unsigned const* candidates,
    unsigned const* gen_offset_of_edges)
{
  unsigned nedges = mesh_count(m, 1);
  unsigned const* tets_of_edges_offsets =
    mesh_ask_up(m, 1, 3)->offsets;
  unsigned const* tets_of_edges =
    mesh_ask_up(m, 1, 3)->adj;
  unsigned const* tets_of_edges_directions =
    mesh_ask_up(m, 1, 3)->directions;
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* verts_of_tets = mesh_ask_down(m, 3, 0);
  double const* coords = mesh_find_tag(m, 0, "coordinates")->d.f64;
  return swap_topology(ent_dim, nedges,
      candidates, gen_offset_of_edges,
      tets_of_edges_offsets, tets_of_edges, tets_of_edges_directions,
      verts_of_edges, verts_of_tets, coords);
}
