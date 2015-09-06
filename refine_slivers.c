#include "refine_slivers.h"
#include "refine_common.h"
#include "mesh.h"
#include "tables.h"
#include "quality.h"
#include "ints.h"
#include <stdlib.h>

unsigned refine_slivers(
    struct mesh** p_m,
    unsigned qual_floor)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned nedges = mesh_count(m, 1);
  unsigned* candidates = malloc(sizeof(unsigned) * nedges);
  unsigned const* elems_of_edges_offsets =
    mesh_ask_up(m, 1, elem_dim)->offsets;
  unsigned const* elems_of_edges =
    mesh_ask_up(m, 1, elem_dim)->adj;
  double* quals = mesh_qualities(m);
  for (unsigned i = 0; i < nedges; ++i) {
    unsigned first_use = elems_of_edges_offsets[i];
    unsigned end_use = elems_of_edges_offsets[i + 1];
    candidates[i] = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned elem = elems_of_edges[j];
      if (quals[elem] < qual_floor)
        candidates[i] = 1;
    }
  }
  free(quals);
  if (!ints_max(candidates, nedges)) {
    free(candidates);
    return 0;
  }
  refine_common(p_m, 1, candidates);
  free(candidates);
  return 1;
}
