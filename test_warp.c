#define _XOPEN_SOURCE 500
#include "mesh.h"
#include "classify_box.h"
#include "refine_by_size.h"
#include "vtk.h"
#include "algebra.h"
#include "warp_to_limit.h"
#include "eval_field.h"
#include <math.h>
#include <stdlib.h>

static double size_fun(double const x[])
{
  return 0.2;
}

static double const max_rot = M_PI / 4.;

static void warp_fun(double const coords[3], double v[])
{
  double x[3];
  double const mid[3] = {.5, .5, 0};
  subtract_vectors(coords, mid, x, 3);
  double polar_a = atan2(x[1], x[0]);
  double polar_r = vector_norm(x, 3);
  double rot_a = max_rot * (2 * (.5 - polar_r));
  if (rot_a < 0)
    rot_a = 0;
  double dest_a = polar_a + rot_a;
  v[0] = cos(dest_a) * polar_r;
  v[1] = sin(dest_a) * polar_r;
  v[2] = 0;
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  while (refine_by_size(&m, size_fun));
  mesh_classify_box(m);
  write_vtk(m, "before.vtu");
  mesh_eval_field(m, "warp", 3, warp_fun);
  mesh_warp_to_limit(m);
  write_vtk(m, "after.vtu");
  free_mesh(m);
}
