#include "swap_qualities.h"
#include "loop.h"     // for malloc
#include "algebra.h"    // for copy_vector
#include "edge_ring.h"  // for find_edge_ring
#include "edge_swap.h"  // for swap_choice, MAX_EDGE_SWAP, choose_edge_swap

void swap_qualities(
    unsigned nedges,
    unsigned* candidates,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    double const* coords,
    double const* elem_quals,
    double** p_qualities,
    unsigned** p_codes,
    unsigned** p_gen_elems_per_edge)
{
  double* out_quals = loop_malloc(sizeof(double) * nedges);
  unsigned* out_codes = loop_malloc(sizeof(unsigned) * nedges);
  unsigned* gen_elems_per_edge = loop_malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    if (!candidates[i])
      continue;
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
      continue;
    }
    double edge_x[2][3];
    copy_vector(coords + edge_v[0] * 3, edge_x[0], 3);
    copy_vector(coords + edge_v[1] * 3, edge_x[1], 3);
    double ring_x[MAX_EDGE_SWAP][3];
    for (unsigned j = 0; j < ring_size; ++j)
      copy_vector(coords + ring_v[j] * 3, ring_x[j], 3);
    struct swap_choice sc = choose_edge_swap(ring_size, edge_x,
        ring_x, old_minq + 1e-10);
    if (sc.quality <= old_minq) {
      candidates[i] = 0;
      gen_elems_per_edge[i] = 0;
    } else {
      out_quals[i] = sc.quality;
      out_codes[i] = sc.code;
      gen_elems_per_edge[i] = 2 * swap_mesh_sizes[ring_size];
    }
  }
  *p_qualities = out_quals;
  *p_codes = out_codes;
  *p_gen_elems_per_edge = gen_elems_per_edge;
}
