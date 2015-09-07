#include "coarsen_slivers.h"
#include "mark_down.h"
#include "quality.h"
#include "ints.h"
#include "mesh.h"
#include "collapse_codes.h"
#include "coarsen_common.h"
#include <stdlib.h>

unsigned coarsen_slivers(
    struct mesh** p_m,
    double quality_floor)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  double* quals = mesh_qualities(m);
  unsigned* slivers = malloc(sizeof(unsigned) * nelems);
  for (unsigned i = 0; i < nelems; ++i)
    slivers[i] = quals[i] < quality_floor;
  unsigned* sliver_offsets = ints_exscan(slivers, nelems);
  free(slivers);
  unsigned* vert_offsets = mesh_mark_down(m, elem_dim, 0,
      sliver_offsets);
  free(sliver_offsets);
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned* col_codes = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    unsigned const* verts_of_edge = verts_of_edges + i * 2;
    col_codes[i] = 0;
    for (unsigned j = 0; j < 2; ++j) {
      unsigned vert = verts_of_edge[j];
      if (vert_offsets[vert] != vert_offsets[vert + 1])
        col_codes[i] |= (1<<j);
    }
  }
  free(vert_offsets);
  if (ints_max(col_codes, nedges) == DONT_COLLAPSE) {
    /* early return #1: no edges are small */
    free(col_codes);
    return 0;
  }
  unsigned ret = coarsen_common(&m, col_codes, 0, 1);
  free(col_codes);
  *p_m = m;
  return ret;
}
