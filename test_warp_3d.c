#include <math.h>            // for atan2, cos, sin, M_PI
#include <stdio.h>           // for printf, fprintf, stderr
#include <stdlib.h>          // for abort
#include "adapt.h"           // for mesh_adapt
#include "algebra.h"         // for subtract_vectors, vector_norm
#include "classify_box.h"    // for mesh_classify_box
#include "eval_field.h"      // for mesh_eval_field
#include "mesh.h"            // for free_mesh, mesh_free_nodal_field, new_bo...
#include "refine_by_size.h"  // for refine_by_size
#include "vtk.h"             // for write_vtk_step, start_vtk_steps
#include "warp_to_limit.h"   // for mesh_warp_to_limit

struct mesh;

static double const warp_qual_floor = 0.2;
static double const good_qual_floor = 0.3;
static double const size_floor = 1. / 3.;
static unsigned const nsliver_layers = 4;
static unsigned const max_ops = 50;

static void size_fun(double const x[], double s[])
{
  (void)x;
  s[0] = 0.2;
}

static double the_rotation = M_PI / 4;

static void warp_fun(double const coords[3], double v[])
{
  double x[3];
  double const mid[3] = {.5, .5, 0};
  subtract_vectors(coords, mid, x, 3);
  x[2] = 0;
  double polar_a = atan2(x[1], x[0]);
  double polar_r = vector_norm(x, 3);
  double rot_a = 0;
  if (polar_r < .5)
    rot_a = the_rotation * (2 * (.5 - polar_r));
  double dest_a = polar_a + rot_a;
  double dst[3];
  dst[0] = cos(dest_a) * polar_r;
  dst[1] = sin(dest_a) * polar_r;
  dst[2] = 0;
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

/*
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
*/

static void warped_adapt(struct mesh** p_m)
{
  static unsigned const n = 6;
  for (unsigned i = 0; i < n; ++i) {
    printf("\n WARP TO LIMIT %u\n", i);
    unsigned done = mesh_warp_to_limit(*p_m, warp_qual_floor);
  //set_size_field(*p_m);
    write_vtk_step(*p_m);
    mesh_adapt(p_m, size_floor, good_qual_floor, nsliver_layers, max_ops);
    if (done)
      return;
  }
  fprintf(stderr, "warped_adapt still not done after %u iters\n", n);
  abort();
}

int main()
{
  struct mesh* m = new_box_mesh(3);
  mesh_eval_field(m, "adapt_size", 1, size_fun);
  while (refine_by_size(&m, 0));
  mesh_classify_box(m);
  start_vtk_steps("warp");
  mesh_eval_field(m, "dye", 1, dye_fun);
  write_vtk_step(m);
  for (unsigned i = 0; i < 2; ++i) {
    printf("\nOUTER DIRECTION %u\n", i);
    for (unsigned j = 0; j < 4; ++j) {
      printf("\nWARP FIELD %u\n", j);
      mesh_eval_field(m, "warp", 3, warp_fun);
      printf("new warp field\n");
      warped_adapt(&m);
      mesh_free_nodal_field(m, "warp");
    }
    the_rotation = -the_rotation;
  }
  free_mesh(m);
}
