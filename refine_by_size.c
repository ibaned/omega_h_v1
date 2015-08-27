#include "refine_by_size.h"
#include "refine_common.h"
#include "measure_edges.h"
#include "ints.h"
#include <stdlib.h>

unsigned refine_by_size(struct mesh** p_m)
{
  struct mesh* m = *p_m;
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned nedges = mesh_count(m, 1);
  double const* sizes = mesh_find_nodal_field(m, "adapt_size")->data;
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, sizes);
  unsigned* candidates = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i)
    candidates[i] = edge_sizes[i] > 1.0;
  free(edge_sizes);
  unsigned something_to_do = ints_max(candidates, nedges);
  if (!something_to_do) {
    free(candidates);
    return 0;
  }
  refine_common(p_m, 1, candidates);
  free(candidates);
  return 1;
}
