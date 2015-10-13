#include "coarsen_by_size.h"

#include "coarsen_common.h"
#include "collapse_codes.h"
#include "loop.h"
#include "measure_edges.h"
#include "mesh.h"
#include "tag.h"

unsigned coarsen_by_size(
    struct mesh** p_m,
    double quality_floor,
    double size_ratio_floor)
{
  struct mesh* m = *p_m;
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  double const* coords = mesh_find_tag(m, 0, "coordinates")->data;
  double const* sizes = mesh_find_tag(m, 0, "adapt_size")->data;
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, sizes);
  unsigned* col_codes = LOOP_MALLOC(unsigned, nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    if (edge_sizes[i] < size_ratio_floor)
      col_codes[i] = COLLAPSE_BOTH;
    else
      col_codes[i] = DONT_COLLAPSE;
  }
  LOOP_FREE(edge_sizes);
  unsigned ret = coarsen_common(&m, col_codes, quality_floor, 0);
  LOOP_FREE(col_codes);
  *p_m = m;
  return ret;
}
