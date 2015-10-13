#include "swap_slivers.h"

#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "swap_common.h"

unsigned swap_slivers(
    struct mesh** p_m,
    double good_qual,
    unsigned nlayers)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned* slivers = mesh_mark_slivers(m, good_qual, nlayers);
  unsigned* candidates = mesh_mark_down(m, elem_dim, 1, slivers);
  LOOP_FREE(slivers);
  unsigned ret = swap_common(p_m, candidates);
  LOOP_FREE(candidates);
  return ret;
}
