#include "refine.hpp"

#include "arrays.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "refine_common.hpp"
#include "size.hpp"

LOOP_KERNEL(refine_candidate,
    double const* edge_sizes,
    unsigned* candidates)
  candidates[i] = edge_sizes[i] > 1.0;
}

unsigned refine_by_size(struct mesh* m, double qual_floor)
{
  double* edge_sizes = mesh_measure_edges_for_adapt(m);
  unsigned nedges = mesh_count(m, 1);
  unsigned* candidates = LOOP_MALLOC(unsigned, nedges);
  LOOP_EXEC(refine_candidate, nedges, edge_sizes, candidates);
  loop_free(edge_sizes);
  mesh_add_tag(m, 1, TAG_U32, "candidate", 1, candidates);
  unsigned ret = refine_common(m, 1, qual_floor, 0);
  return ret;
}

void uniformly_refine(struct mesh* m)
{
  unsigned nedges = mesh_count(m, 1);
  mesh_add_tag(m, 1, TAG_U32, "candidate", 1, uints_filled(nedges, 1));
  refine_common(m, 1, 0.0, 0);
}
