#include "refine_qualities.h"
#include "quality.h"
#include "average.h"
#include "tables.h"
#include <assert.h>
#include <stdlib.h>

double* refine_qualities(
    unsigned elem_dim,
    unsigned split_dim,
    unsigned nsplit,
    unsigned const* split_verts,
    unsigned const* split_elem_offsets,
    unsigned const* split_elem_edges,
    unsigned const* split_elem_directions,
    unsigned const* elem_verts,
    double const* coords)
{
  assert(elem_dim <= 3);
  assert(elem_dim >= split_dim);
  unsigned nelem_vert = the_down_degrees[elem_dim][0];
  unsigned base_dim = elem_dim - 1;
  unsigned nbase_vert = the_down_degrees[base_dim][0];
  unsigned opposite_dim = get_opposite_dim(elem_dim, base_dim);
  assert(opposite_dim == 0);
  unsigned nsplit_vert = the_down_degrees[split_dim][0];
  assert(nsplit_vert);
  unsigned const* const* any_split_opposite =
    the_canonical_orders[elem_dim][split_dim][opposite_dim];
  unsigned const* const* any_base_vert =
    the_canonical_orders[elem_dim][base_dim][0];
  unsigned const* opposite_base = the_opposite_orders[elem_dim][opposite_dim];
  quality_function qf = the_quality_functions[elem_dim];
  double* out = malloc(sizeof(double) * nsplit);
  for (unsigned i = 0; i < nsplit; ++i) {
    unsigned first_elem = split_elem_offsets[i];
    unsigned last_elem = split_elem_offsets[i + 1];
    if (last_elem == first_elem)
      continue;
    double elem_coords[MAX_DOWN][3];
    average_element_field(
        nsplit_vert,
        split_verts + i * nsplit_vert,
        coords,
        3,
        elem_coords[nelem_vert - 1]);
    double minq = 1;
    for (unsigned j = first_elem; j < last_elem; ++j) {
      unsigned elem = split_elem_edges[j];
      unsigned const* elem_vert = elem_verts + elem * nelem_vert;
      unsigned direction = split_elem_directions[j];
      unsigned const* split_opposite = any_split_opposite[direction];
      for (unsigned k = 0; k < nsplit_vert; ++k) {
        unsigned opposite = split_opposite[k];
        unsigned base = opposite_base[opposite];
        unsigned const* base_vert = any_base_vert[base];
        for (unsigned l = 0; l < nbase_vert; ++l) {
          unsigned vert = elem_vert[base_vert[k]];
          for (unsigned m = 0; m < 3; ++m)
            elem_coords[l][m] = coords[vert * 3 + m];
        }
        double q = qf(elem_coords);
        if (q < minq)
          minq = q;
      }
    }
    out[i] = minq;
  }
  return out;
}
