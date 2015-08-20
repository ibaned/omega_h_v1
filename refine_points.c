#include "refine_points.h"
#include "inside.h"
#include "tables.h"
#include "algebra.h"
#include <stdlib.h>
#include <assert.h>

static unsigned transfer_elem_points(
    unsigned elem_dim,
    unsigned const* verts_of_new_elem,
    double const* coords_of_verts,
    unsigned npts_of_src_elem,
    unsigned const* pts_of_src_elem,
    double const* coords_of_pts,
    unsigned pt_buf[])
{
  unsigned nverts_per_elem = the_down_degrees[elem_dim][0];
  double coords_of_elem[4][3];
  for (unsigned i = 0; i < nverts_per_elem; ++i) {
    unsigned vert = verts_of_new_elem[i];
    copy_vector(coords_of_verts + vert * 3, coords_of_elem[i], 3);
  }
  unsigned nin = 0;
  for (unsigned i = 0; i < npts_of_src_elem; ++i) {
    unsigned pt = pts_of_src_elem[i];
    double b[4];
    the_inside_functions[elem_dim](coords_of_elem, coords_of_pts + pt * 3, b);
    if (is_in_simplex(elem_dim, b))
      pt_buf[nin++] = pt;
  }
  return nin;
}

void refine_points(
    unsigned elem_dim,
    unsigned nsrc_elems,
    unsigned nnew_elems,
    unsigned const* verts_of_new_elems,
    double const* coords_of_verts,
    unsigned const* pts_of_src_elems_offsets,
    unsigned const* pts_of_src_elems,
    double const* coords_of_pts,
    unsigned** p_pts_of_new_elems_offsets,
    unsigned** p_pts_of_new_elems)
{
  assert(nnew_elems % nsrc_elems == 0);
  unsigned new_elems_per_src_elem = nnew_elems / nsrc_elems;
  unsigned npts_of_src_elems = pts_of_src_elems_offsets[nsrc_elems];
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* pts_of_new_elems = malloc(sizeof(unsigned) * npts_of_src_elems);
  unsigned* pts_of_new_elems_offsets = malloc(sizeof(unsigned) * (nnew_elems + 1));
  pts_of_new_elems_offsets[0] = 0;
  for (unsigned i = 0; i < nsrc_elems; ++i) {
    unsigned npts_of_src_elem = pts_of_src_elems_offsets[i + 1] -
      pts_of_src_elems_offsets[i];
    if (!npts_of_src_elem)
      continue;
    unsigned const* pts_of_src_elem  = pts_of_src_elems +
      pts_of_src_elems_offsets[i];
    unsigned* buf = pts_of_new_elems + pts_of_src_elems_offsets[i];
    for (unsigned j = 0; j < nnew_elems; ++j) {
      unsigned new_elem = i * new_elems_per_src_elem + j;
      unsigned const* verts_of_new_elem = verts_of_new_elems +
        new_elem * verts_per_elem;
      unsigned degree = transfer_elem_points(elem_dim, verts_of_new_elem,
          coords_of_verts, npts_of_src_elem, pts_of_src_elem, coords_of_pts,
          buf);
      buf += degree;
      pts_of_new_elems_offsets[new_elem + 1] = buf - pts_of_new_elems;
    }
    assert(buf == pts_of_new_elems + pts_of_src_elems_offsets[i + 1]);
  }
  *p_pts_of_new_elems_offsets = pts_of_new_elems_offsets;
  *p_pts_of_new_elems = pts_of_new_elems;
}
