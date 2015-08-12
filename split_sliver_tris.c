#include "split_sliver_tris.h"
#include "bad_elem_keys.h"
#include "derive_edges.h"
#include "ints.h"
#include <assert.h>
#include <stdlib.h>

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
  unsigned nelems = *p_nelems;
  unsigned nverts = *p_nverts;
  unsigned const* verts_of_elems = *p_verts_of_elems;
  unsigned ntris = *p_nelems;
  unsigned const* verts_of_tris = *p_verts_of_elems;
  double const* coords = *p_coords;
  unsigned* bad_tris;
  unsigned* key_of_tris;
  bad_elem_keys(2, ntris, verts_of_tris, coords,
      SLIVER_ELEM, qual_floor, edge_ratio_floor,
      &bad_tris, &key_of_tris);
  unsigned something_to_do = ints_max(bad_tris, ntris);
  if (!something_to_do) {
    free(bad_tris);
    free(key_of_tris);
    return 0;
  }
  unsigned nedges;
  unsigned* verts_of_edges;
  derive_edges(elem_dim, nelems, nverts, verts_of_elems,
      &nedges, &verts_of_edges);
  return 1;
}
