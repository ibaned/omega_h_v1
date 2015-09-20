#include "coarsen_by_size.h"
#include "coarsen_common.h"
#include "mesh.h"  // for mesh_ask_down, mesh_count
#include "field.h"
#include "measure_edges.h"
#include "collapse_codes.h"
#include "loop.h"

unsigned coarsen_by_size(
    struct mesh** p_m,
    double quality_floor,
    double size_ratio_floor)
{
  struct mesh* m = *p_m;
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  double const* coords = mesh_find_field(m, 0, "coordinates")->data;
  double const* sizes = mesh_find_field(m, 0, "adapt_size")->data;
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, sizes);
  unsigned* col_codes = loop_malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    if (edge_sizes[i] < size_ratio_floor)
      col_codes[i] = COLLAPSE_BOTH;
    else
      col_codes[i] = DONT_COLLAPSE;
  }
  loop_free(edge_sizes);
  unsigned ret = coarsen_common(&m, col_codes, quality_floor, 0);
  loop_free(col_codes);
  *p_m = m;
  return ret;
}
