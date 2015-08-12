#include "split_sliver_tris.h"
#include <assert.h>

int split_sliver_tris(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    double qual_floor,
    double edge_ratio_floor)
{
  /* TODO: accept a tet mesh */
  assert(elem_dim == 2);
  unsigned ntris = *p_nelems;
  unsigned const* verts_of_tris = *p_verts_of_elems;
  (void) ntris;
  (void) verts_of_tris;
  return 1;
}
