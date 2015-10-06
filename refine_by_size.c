#include "refine_by_size.h"

#include "loop.h"         // for free, malloc
#include "measure_edges.h"  // for measure_edges
#include "mesh.h"           // for mesh_find_field, mesh_ask_down, mes...
#include "refine_common.h"  // for refine_common
#include "tag.h"

unsigned refine_by_size(struct mesh** p_m, double qual_floor)
{
  struct mesh* m = *p_m;
  double const* coords = mesh_find_tag(m, 0, "coordinates")->data;
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned nedges = mesh_count(m, 1);
  double const* sizes = mesh_find_tag(m, 0, "adapt_size")->data;
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, sizes);
  unsigned* candidates = loop_malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i)
    candidates[i] = edge_sizes[i] > 1.0;
  loop_free(edge_sizes);
  unsigned ret = refine_common(p_m, 1, candidates, qual_floor, 0);
  loop_free(candidates);
  return ret;
}
