#include "form_cloud.h"
#include "mesh.h"
#include "cloud.h"
#include "loop.h"
#include "size.h"
#include "doubles.h"
#include <math.h>
#include <assert.h>

#include <stdio.h>

static void form_elem_cloud_2d(
    unsigned elem,
    unsigned const* elem_pt_offsets,
    unsigned const* verts_of_elems,
    double const* vert_coords,
    double* coords_of_pts,
    unsigned* elem_of_pts)
{
  fprintf(stderr, "forming cloud for element %u\n", elem);
  unsigned elem_pts = (elem_pt_offsets[elem + 1]
                     - elem_pt_offsets[elem]);
  fprintf(stderr, "%u points in this element\n", elem_pts);
  if (!elem_pts)
    return;
  unsigned res;
  for (res = 1; ((res * (res + 1)) / 2) < elem_pts; ++res);
  fprintf(stderr, "resolution %u\n", res);
  double dxi = 1.0 / (res + 2);
  fprintf(stderr, "dxi %f\n", dxi);
  unsigned elem_pt = 0;
  unsigned pt_offset = elem_pt_offsets[elem];
  double const* elem_x[3];
  unsigned const* verts_of_elem = verts_of_elems + elem * 3;
  for (unsigned i = 0; i < 3; ++i)
    elem_x[i] = vert_coords + verts_of_elem[i] * 3;
  for (unsigned i = 0; i < res; ++i)
  for (unsigned j = 0; j < (res - i); ++j) {
    if (elem_pt == elem_pts)
      return;
    fprintf(stderr, "local point %u\n", elem_pt);
    double xi[3];
    xi[0] = (i + 1) * dxi;
    xi[1] = (j + 1) * dxi;
    xi[2] = 1.0 - xi[0] - xi[1];
    fprintf(stderr, "parent coordinates %f %f %f\n", xi[0], xi[1], xi[2]);
    double* coords_of_pt = coords_of_pts + (pt_offset + elem_pt) * 3;
    for (unsigned k = 0; k < 3; ++k)
      coords_of_pt[k] = 0;
    for (unsigned k = 0; k < 3; ++k)
    for (unsigned l = 0; l < 3; ++l)
      coords_of_pt[l] += elem_x[k][l] * xi[k];
    elem_of_pts[pt_offset + elem_pt] = elem;
    ++elem_pt;
  }
}

struct cloud* form_cloud(struct mesh* m)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  double const* elem_density = mesh_find_tag(m, dim, "cloud_density")->data;
  double const* elem_size = mesh_element_sizes(m);
  double* pts_of_elems = loop_malloc(sizeof(double) * nelems);
  for (unsigned i = 0; i < nelems; ++i) {
    pts_of_elems[i] = elem_size[i] * elem_density[i];
    fprintf(stderr, "%f real points for element %u\n", pts_of_elems[i], i);
  }
  double* real_offset = doubles_exscan(pts_of_elems, nelems);
  loop_free(pts_of_elems);
  unsigned* uint_offset = loop_malloc(sizeof(unsigned) * (nelems + 1));
  for (unsigned i = 0; i <= nelems; ++i) {
    uint_offset[i] = (unsigned) (floor(real_offset[i]));
  }
  loop_free(real_offset);
  unsigned npts = uint_offset[nelems];
  fprintf(stderr, "%u total cloud points\n", npts);
  struct cloud* c = new_cloud(npts);
  unsigned* pt_elem = loop_malloc(sizeof(unsigned) * npts);
  double* pt_coords = loop_malloc(sizeof(double) * npts * 3);
  unsigned const* verts_of_elems = mesh_ask_down(m, dim, 0);
  double const* vert_coords = mesh_find_tag(m, 0, "coordinates")->data;
  assert(dim == 2);
  for (unsigned i = 0; i < nelems; ++i)
    form_elem_cloud_2d(i, uint_offset, verts_of_elems, vert_coords,
        pt_coords, pt_elem);
  cloud_add_tag(c, TAG_F64, "coordinates", 3, pt_coords);
  cloud_add_tag(c, TAG_U32, "mesh_elem", 1, pt_elem);
  return c;
}
