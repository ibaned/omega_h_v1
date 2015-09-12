#include "coarsen_by_size.h"
#include <stdlib.h>          // for free, malloc
#include "coarsen_common.h"  // for coarsen_common
#include "collapse_codes.h"  // for ::COLLAPSE_BOTH, ::DONT_COLLAPSE
#include "field.h"           // for const_field
#include "measure_edges.h"   // for measure_edges
#include "mesh.h"            // for mesh_find_nodal_field, mesh_ask_down

unsigned coarsen_by_size(
    struct mesh** p_m,
    double quality_floor,
    double size_ratio_floor,
    unsigned require_better)
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
  unsigned ret = coarsen_common(&m, col_codes, quality_floor, require_better);
  free(col_codes);
  *p_m = m;
  return ret;
}
