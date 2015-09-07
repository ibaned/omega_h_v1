#include "refine_slivers.h"
#include "refine_common.h"
#include "mesh.h"
#include "tables.h"
#include "quality.h"
#include "ints.h"
#include "mark_down.h"
#include <stdlib.h>

unsigned refine_slivers(
    struct mesh** p_m,
    double good_qual,
    double valid_qual,
    unsigned require_better)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned* marked_elems = malloc(sizeof(unsigned) * nelems);
  double* quals = mesh_qualities(m);
  for (unsigned i = 0; i < nelems; ++i)
    marked_elems[i] = quals[i] < good_qual;
  free(quals);
  unsigned* candidates = mesh_mark_down(m, elem_dim, 1, marked_elems);
  free(marked_elems);
  unsigned ret = refine_common(p_m, 1, candidates, valid_qual, require_better);
  free(candidates);
  return ret;
}
