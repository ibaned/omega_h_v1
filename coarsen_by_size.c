#include "coarsen_by_size.h"
#include "coarsen_common.h"
#include "measure_edges.h"
#include "collapse_codes.h"
#include "ints.h"
#include <stdlib.h>

unsigned coarsen_by_size(
    struct mesh** p_m,
    double quality_floor,
    double size_ratio_floor)
{
  struct mesh* m = *p_m;
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  double const* sizes = mesh_find_nodal_field(m, "adapt_size")->data;
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, sizes);
  unsigned* col_codes = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    if (edge_sizes[i] < size_ratio_floor)
      col_codes[i] = COLLAPSE_BOTH;
    else
      col_codes[i] = DONT_COLLAPSE;
  }
  free(edge_sizes);
  if (ints_max(col_codes, nedges) == DONT_COLLAPSE) {
    /* early return #1: no edges are small */
    free(col_codes);
    return 0;
  }
  unsigned ret = coarsen_common(&m, col_codes, quality_floor, size_ratio_floor, 0);
  free(col_codes);
  *p_m = m;
  return ret;
}
