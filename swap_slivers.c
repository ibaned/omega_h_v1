#include "swap_slivers.h"
#include "swap_common.h"
#include "mesh.h"
#include "mark.h"
#include "loop.h"

unsigned swap_slivers(
    struct mesh** p_m,
    double good_qual,
    double valid_qual,
    unsigned nlayers)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned* slivers = mesh_mark_slivers(m, good_qual, nlayers);
  unsigned* candidates = mesh_mark_down(m, elem_dim, 1, slivers);
  loop_free(slivers);
  unsigned ret = swap_common(p_m, candidates, 1.0, valid_qual, 1);
  loop_free(candidates);
  return ret;
}
