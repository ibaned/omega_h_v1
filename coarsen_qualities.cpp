#include "coarsen_qualities.hpp"

#include "algebra.hpp"
#include "arrays.hpp"
#include "collapse_codes.hpp"
#include "loop.hpp"
#include "quality.hpp"
#include "tables.hpp"

LOOP_KERNEL(coarsen_quality,
    unsigned elem_dim,
    unsigned base_dim,
    unsigned const* verts_of_edges,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    unsigned const* elems_of_verts_directions,
    unsigned const* verts_of_elems,
    unsigned verts_per_elem,
    double const* coords,
    double const* elem_quals,
    double quality_floor,
    unsigned* col_codes,
    double* out)

  if (col_codes[i] == DONT_COLLAPSE)
    return;
  unsigned const* const* elem_verts_of_bases =
      the_canonical_orders[elem_dim][base_dim][0];
  unsigned const* elem_bases_opp_verts =
      the_opposite_orders[elem_dim][0];
  unsigned const* verts_of_edge = verts_of_edges + i * 2;
  unsigned require_better = (elem_quals != 0);
  for (unsigned j = 0; j < 2; ++j) {
    if (!collapses(col_codes[i], j))
      continue;
    double minq = 1;
    double old_minq = 1;
    unsigned edge_direction = j;
    unsigned col_vert = verts_of_edge[edge_direction];
    unsigned gen_vert = verts_of_edge[1 - edge_direction];
    double elem_x[MAX_DOWN][3];
    copy_vector(coords + gen_vert * 3, elem_x[verts_per_elem - 1], 3);
    unsigned first_use = elems_of_verts_offsets[col_vert];
    unsigned end_use = elems_of_verts_offsets[col_vert + 1];
    for (unsigned k = first_use; k < end_use; ++k) {
      unsigned elem = elems_of_verts[k];
      unsigned const* verts_of_elem = verts_of_elems + elem * verts_per_elem;
      unsigned is_adjacent_to_edge = 0;
      if (require_better) {
        double old_q = elem_quals[elem];
        if (old_q < old_minq)
          old_minq = old_q;
      }
      for (unsigned l = 0; l < verts_per_elem; ++l)
        if (verts_of_elem[l] == gen_vert)
          is_adjacent_to_edge = 1;
      if (is_adjacent_to_edge)
        continue;
      unsigned elem_direction = elems_of_verts_directions[k];
      unsigned base = elem_bases_opp_verts[elem_direction];
      unsigned const* elem_verts_of_base = elem_verts_of_bases[base];
      for (unsigned l = 0; l < (verts_per_elem - 1); ++l) {
        unsigned vert = verts_of_elem[elem_verts_of_base[l]];
        copy_vector(coords + vert * 3, elem_x[l], 3);
      }
      double q = element_quality(elem_dim, elem_x);
      if (q < minq)
        minq = q;
    }
    if ((minq < quality_floor) ||
        (require_better && (minq <= old_minq))) {
      col_codes[i] = dont_collapse(col_codes[i], j);
    } else {
      out[i * 2 + j] = minq;
    }
  }
}

double* coarsen_qualities(
    unsigned elem_dim,
    unsigned nedges,
    unsigned* col_codes,
    unsigned const* verts_of_elems,
    unsigned const* verts_of_edges,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    unsigned const* elems_of_verts_directions,
    double const* coords,
    double quality_floor,
    double const* elem_quals)
{
  if (elem_dim < 2)
    return filled_array(nedges * 2, 1.0);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned base_dim = elem_dim - 1;
  double* out = LOOP_MALLOC(double, nedges * 2);
  LOOP_EXEC(coarsen_quality, nedges,
      elem_dim,
      base_dim,
      verts_of_edges,
      elems_of_verts_offsets,
      elems_of_verts,
      elems_of_verts_directions,
      verts_of_elems,
      verts_per_elem,
      coords,
      elem_quals,
      quality_floor,
      col_codes,
      out);
  return out;
}
