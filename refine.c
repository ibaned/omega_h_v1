#include "refine.h"

#include "arrays.h"
#include "loop.h"
#include "mesh.h"
#include "refine_common.h"
#include "size.h"

unsigned refine_by_size(struct mesh** p_m, double qual_floor)
{
  struct mesh* m = *p_m;
  double* edge_sizes = mesh_measure_edges_for_adapt(m);
  unsigned nedges = mesh_count(m, 1);
  unsigned* candidates = LOOP_MALLOC(unsigned, nedges);
  for (unsigned i = 0; i < nedges; ++i)
    candidates[i] = edge_sizes[i] > 1.0;
  loop_free(edge_sizes);
  mesh_add_tag(m, 1, TAG_U32, "candidate", 1, candidates);
  unsigned ret = refine_common(p_m, 1, qual_floor, 0);
  return ret;
}

void uniformly_refine(struct mesh** p_m)
{
  struct mesh* m = *p_m;
  unsigned nedges = mesh_count(m, 1);
  mesh_add_tag(m, 1, TAG_U32, "candidate", 1, uints_filled(nedges, 1));
  refine_common(p_m, 1, 0.0, 0);
}
