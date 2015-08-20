#include "embed_points.h"
#include "inside.h"
#include "algebra.h"
#include "tables.h"
#include "ints.h"
#include <stdlib.h>

static unsigned embed_elem_points(
    unsigned elem_dim,
    inside_function insf,
    unsigned verts_per_elem,
    unsigned npts,
    unsigned const* verts_of_elem,
    double const* coords_of_verts,
    double const* coords_of_pts,
    unsigned pt_buf[])
{
  double coords_of_elem[4][3];
  for (unsigned i = 0; i < verts_per_elem; ++i) {
    unsigned vert = verts_of_elem[i];
    copy_vector(coords_of_verts + vert * 3, coords_of_elem[i], 3);
  }
  unsigned nin = 0;
  for (unsigned j = 0; j < npts; ++j) {
    double b[4];
    insf(coords_of_elem, coords_of_pts + j * 3, b);
    if (is_in_simplex(elem_dim, b)) {
      if (pt_buf)
        pt_buf[nin] = j;
      ++nin;
    }
  }
  return nin;
}

void embed_points(
    unsigned elem_dim,
    unsigned nelems,
    unsigned npts,
    unsigned const* verts_of_elems,
    double const* coords_of_verts,
    double const* coords_of_pts,
    unsigned** p_pts_of_elems_offsets,
    unsigned** p_pts_of_elems)
{
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* degrees = malloc(sizeof(unsigned) * nelems);
  inside_function insf = the_inside_functions[elem_dim];
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    degrees[i] = embed_elem_points(elem_dim, insf, verts_per_elem, npts,
        verts_of_elem, coords_of_verts, coords_of_pts, 0);
  }
  unsigned* pts_of_elems_offsets = ints_exscan(degrees, nelems);
  free(degrees);
  unsigned* pts_of_elems = malloc(sizeof(unsigned) * npts);
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    unsigned* pts_of_elem = pts_of_elems + pts_of_elems_offsets[i];
    degrees[i] = embed_elem_points(elem_dim, insf, verts_per_elem, npts,
        verts_of_elem, coords_of_verts, coords_of_pts, pts_of_elem);
  }
  *p_pts_of_elems_offsets = pts_of_elems_offsets;
  *p_pts_of_elems = pts_of_elems;
}
