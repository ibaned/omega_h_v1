#define _XOPEN_SOURCE 500
#include "mesh.h"
#include "classify_box.h"
#include "refine_by_size.h"
#include "vtk.h"
#include "algebra.h"
#include "warp_to_limit.h"
#include "eval_field.h"
#include "split_sliver_tris.h"
#include "coarsen_by_size.h"
#include "quality.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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
  double dst[3];
  dst[0] = cos(dest_a) * polar_r;
  dst[1] = sin(dest_a) * polar_r;
  dst[2] = 0;
  subtract_vectors(dst, x, v, 3);
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  while (refine_by_size(&m, size_fun));
  mesh_classify_box(m);
  start_vtk_steps("warp");
  write_vtk_step(m);
  while (1) {
    mesh_eval_field(m, "warp", 3, warp_fun);
    unsigned ok = mesh_warp_to_limit(m, 0.1);
    mesh_free_nodal_field(m, "warp");
    write_vtk_step(m);
    if (!ok)
      break;
  }
  while (1) {
    unsigned did_stuff = 0;
    did_stuff |= coarsen_by_size(&m, size_fun,
        mesh_min_quality(m), 1.0 / 3.0);
    write_vtk_step(m);
    did_stuff |= refine_by_size(&m, size_fun);
    write_vtk_step(m);
    did_stuff |= split_sliver_tris(&m, 0.4, 1.0 / 5.0);
    write_vtk_step(m);
    if (!did_stuff)
      break;
  }
  free_mesh(m);
}
