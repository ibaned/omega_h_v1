#include "refine_slivers.h"
#include "refine_common.h"
#include "mesh.h"
#include "tables.h"
#include "quality.h"
#include "ints.h"
#include "mark.h"
#include <stdlib.h>

unsigned refine_slivers(
    struct mesh** p_m,
    unsigned src_dim,
    double good_qual,
    double valid_qual,
    unsigned require_better,
    unsigned nlayers)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned* slivers = mesh_mark_slivers(m, good_qual, nlayers);
  unsigned* candidates = mesh_mark_down(m, elem_dim, src_dim, slivers);
  free(slivers);
  unsigned ret = refine_common(p_m, src_dim, candidates, valid_qual, require_better);
  free(candidates);
  return ret;
}
