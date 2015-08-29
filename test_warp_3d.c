#define _XOPEN_SOURCE 500
#include "mesh.h"
#include "classify_box.h"
#include "refine_by_size.h"
#include "vtk.h"
#include "algebra.h"
#include "warp_to_limit.h"
#include "eval_field.h"
#include "split_slivers.h"
#include "coarsen_by_size.h"
#include "quality.h"
#include "element_gradients.h"
#include "recover_by_volume.h"
#include "size_from_hessian.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static void size_fun(double const x[], double s[])
{
  s[0] = 0.1;
}

static double the_rotation = M_PI / 4.;

static void warp_fun(double const coords[3], double v[])
{
  double x[3];
  double const mid[3] = {.5, .5, 0};
  subtract_vectors(coords, mid, x, 3);
  double polar_a = atan2(x[1], x[0]);
  double polar_r = vector_norm(x, 3);
  double rot_a = 0;
  if (polar_r < .5)
    rot_a = the_rotation * (2 * (.5 - polar_r));
  double dest_a = polar_a + rot_a;
  double dst[3];
  dst[0] = cos(dest_a) * polar_r;
  dst[1] = sin(dest_a) * polar_r;
  dst[2] = x[2];
  subtract_vectors(dst, x, v, 3);
}

static void dye_fun(double const coords[3], double v[])
{
  double x[3];
  double const l[3] = {.25, .5, .5};
  double const r[3] = {.75, .5, .5};
  double dir = 1;
  subtract_vectors(coords, l, x, 3);
  if (vector_norm(x, 3) > .25) {
    dir = -1;
    subtract_vectors(coords, r, x, 3);
  }
  if (vector_norm(x, 3) > .25) {
    v[0] = 0;
    return;
  }
  v[0] = 4 * dir * (.25 - vector_norm(x, 3));
}

static void set_size_field(struct mesh* m)
{
  mesh_free_nodal_field(m, "adapt_size");
  mesh_element_gradients(m, "dye");
  mesh_recover_by_volume(m, "grad_dye");
  mesh_free_elem_field(m,  "grad_dye");
  mesh_element_gradients(m, "rcov_grad_dye");
  mesh_free_nodal_field(m,  "rcov_grad_dye");
  mesh_recover_by_volume(m, "grad_rcov_grad_dye");
  mesh_free_elem_field(m,  "grad_rcov_grad_dye");
  double weight = 0.05 / 100.0;
  mesh_size_from_hessian(m, "rcov_grad_rcov_grad_dye", &weight, 0.05, 0.1);
  mesh_free_nodal_field(m, "rcov_grad_rcov_grad_dye");
}

static void adapt(struct mesh** p_m)
{
  struct mesh* m = *p_m;
  set_size_field(m);
  write_vtk_step(m);
  while (refine_by_size(&m)) {
    printf("refine\n");
    write_vtk_step(m);
  }
  while (coarsen_by_size(&m, mesh_min_quality(m), 1.0 / 3.0)) {
    printf("coarsen\n");
    write_vtk_step(m);
  }
  while (split_slivers(2, &m, 0.4, 1.0 / 5.0)) {
    printf("sliver\n");
    write_vtk_step(m);
    coarsen_by_size(&m, mesh_min_quality(m), 1.0 / 3.0);
    printf("sliver coarsen\n");
    write_vtk_step(m);
  }
  *p_m = m;
}

static void warped_adapt(struct mesh** p_m)
{
  unsigned done = 0;
  do {
    done = mesh_warp_to_limit(*p_m, 0.1);
    adapt(p_m);
  } while (!done);
}

int main()
{
  struct mesh* m = new_box_mesh(3);
  mesh_eval_field(m, "adapt_size", 1, size_fun);
  while (refine_by_size(&m));
  mesh_classify_box(m);
  start_vtk_steps("warp");
  mesh_eval_field(m, "dye", 1, dye_fun);
  adapt(&m);
//for (unsigned i = 0; i < 0; ++i) {
//  for (unsigned j = 0; j < 0; ++j) {
//    mesh_eval_field(m, "warp", 3, warp_fun);
//    printf("new warp field\n");
//    warped_adapt(&m);
//    mesh_free_nodal_field(m, "warp");
//  }
//  the_rotation = -the_rotation;
//}
  free_mesh(m);
}
