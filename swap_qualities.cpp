#include "swap_qualities.hpp"

#include "algebra.hpp"
#include "edge_ring.hpp"
#include "edge_swap.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "quality.hpp"

LOOP_KERNEL(swap_quality,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    double const* coords,
    double const* elem_quals,
    unsigned const* owned_edges,
    unsigned* candidates,
    double* out_quals,
    unsigned* ring_sizes)

  if (!candidates[i])
    return;
  if (owned_edges && !owned_edges[i])
    return;
  unsigned first_use = tets_of_edges_offsets[i];
  unsigned end_use = tets_of_edges_offsets[i + 1];
  double old_minq = 1;
  for (unsigned j = first_use; j < end_use; ++j) {
    unsigned tet = tets_of_edges[j];
    double tet_q = elem_quals[tet];
    if (tet_q < old_minq)
      old_minq = tet_q;
  }
  unsigned edge_v[2];
  unsigned ring_v[MAX_EDGE_SWAP];
  unsigned ring_size = find_edge_ring(i,
      tets_of_edges_offsets, tets_of_edges, tets_of_edges_directions,
      verts_of_edges, verts_of_tets, edge_v, ring_v);
  if (ring_size > MAX_EDGE_SWAP) {
    candidates[i] = 0;
    return;
  }
  double edge_x[2][3];
  copy_vector(coords + edge_v[0] * 3, edge_x[0], 3);
  copy_vector(coords + edge_v[1] * 3, edge_x[1], 3);
  double ring_x[MAX_EDGE_SWAP][3];
  for (unsigned j = 0; j < ring_size; ++j)
    copy_vector(coords + ring_v[j] * 3, ring_x[j], 3);
  struct swap_choice sc = choose_edge_swap(ring_size, edge_x, ring_x);
  if (sc.quality > old_minq) {
    out_quals[i] = sc.quality;
    ring_sizes[i] = ring_size;
  } else {
    candidates[i] = 0;
  }
}

static void swap_qualities(
    unsigned nedges,
    unsigned* candidates,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    double const* coords,
    double const* elem_quals,
    unsigned const* owned_edges,
    double** p_qualities,
    unsigned** p_ring_sizes)
{
  double* out_quals = LOOP_MALLOC(double, nedges);
  unsigned* ring_sizes = LOOP_MALLOC(unsigned, nedges);
  LOOP_EXEC(swap_quality, nedges,
      tets_of_edges_offsets,
      tets_of_edges,
      tets_of_edges_directions,
      verts_of_edges,
      verts_of_tets,
      coords,
      elem_quals,
      owned_edges,
      candidates,
      out_quals,
      ring_sizes);
  *p_qualities = out_quals;
  *p_ring_sizes = ring_sizes;
}

void mesh_swap_qualities(
    struct mesh* m,
    unsigned* candidates,
    double** p_qualities,
    unsigned** p_ring_sizes)
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
  double* elem_quals = mesh_qualities(m);
  unsigned* owned_edges = 0;
  if (mesh_is_parallel(m))
    owned_edges = mesh_get_owned(m, 1);
  swap_qualities(nedges, candidates,
      tets_of_edges_offsets, tets_of_edges, tets_of_edges_directions,
      verts_of_edges, verts_of_tets,
      coords, elem_quals, owned_edges,
      p_qualities, p_ring_sizes);
  loop_free(elem_quals);
  loop_free(owned_edges);
}
