#include "refine_by_size.h"
#include <stdlib.h>         // for free, malloc
#include "field.h"          // for const_field
#include "measure_edges.h"  // for measure_edges
#include "mesh.h"           // for mesh_find_nodal_field, mesh_ask_down, mes...
#include "refine_common.h"  // for refine_common

unsigned refine_by_size(struct mesh** p_m, double qual_floor)
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
  unsigned ret = refine_common(p_m, 1, candidates, qual_floor, 0);
  free(candidates);
  return ret;
}
