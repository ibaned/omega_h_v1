#define _XOPEN_SOURCE 500
#include "mesh.h"
#include "classify_box.h"
#include "refine_by_size.h"
#include "vtk.h"
#include "algebra.h"
#include "warp_to_limit.h"
#include <math.h>
#include <stdlib.h>

static double size_fun(double const x[])
{
  return 0.2;
}

static double* set_warps(struct mesh* m, double max_rot)
{
  unsigned nverts = mesh_count(m, 0);
  double* warps = malloc(sizeof(double) * 3 * nverts);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  double const mid[3] = {.5, .5, 0};
  for (unsigned i = 0; i < nverts; ++i) {
    double x[3];
    subtract_vectors(coords + i * 3, mid, x, 3);
    double polar_a = atan2(x[1], x[0]);
    double polar_r = vector_norm(x, 3);
    double rot_a = max_rot * (2 * (.5 - polar_r));
    if (rot_a < 0)
      rot_a = 0;
    double* v = warps + i * 3;
    double dest_a = polar_a + rot_a;
    v[0] = cos(dest_a) * polar_r;
    v[1] = sin(dest_a) * polar_r;
    v[2] = 0;
  }
  return warps;
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  while (refine_by_size(&m, size_fun));
  mesh_add_nodal_label(m, "class_dim",
      classify_box(mesh_dim(m), mesh_count(m, 0),
        mesh_find_nodal_field(m, "coordinates")->data));
  write_vtk(m, "before.vtu");
  double* warps = set_warps(m, 2 * M_PI);
  double* new_coords;
  warp_to_limit(2, mesh_count(m, 2), mesh_count(m, 0),
      mesh_ask_down(m, 2, 0), mesh_find_nodal_field(m, "coordinates")->data,
      warps, &new_coords, 0);
  free(warps);
  mesh_free_nodal_field(m, "coordinates");
  mesh_add_nodal_field(m, "coordinates", 3, new_coords);
  write_vtk(m, "after.vtu");
  free_mesh(m);
}
