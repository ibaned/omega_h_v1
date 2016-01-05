#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "adapt.h"
#include "algebra.h"
#include "classify_box.h"
#include "derive_model.h"
#include "element_gradients.h"
#include "eval_field.h"
#include "mesh.h"
#include "recover_by_volume.h"
#include "refine_by_size.h"
#include "size_from_hessian.h"
#include "vtk.h"
#include "warp_to_limit.h"

struct mesh;

static double const warp_qual_floor = 0.1;
static double const good_qual_floor = 0.4;
static double const size_floor = 1. / 3.;
static unsigned const nsliver_layers = 0;
static unsigned const max_ops = 50;

static void size_fun(double const* x, double* s)
{
  (void)x;
  s[0] = 0.1;
}

static double the_rotation = PI / 4.;

static void warp_fun(double const* coords, double* v)
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
  dst[2] = 0;
  subtract_vectors(dst, x, v, 3);
}

static void dye_fun(double const coords[3], double* v)
{
  double x[3];
  double const l[3] = {.25, .5, 0};
  double const r[3] = {.75, .5, 0};
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
  mesh_free_tag(m, 0, "adapt_size");
  mesh_element_gradients(m, "dye");
  mesh_recover_by_volume(m, "grad_dye");
  mesh_free_tag(m, mesh_dim(m),  "grad_dye");
  mesh_element_gradients(m, "grad_dye");
  mesh_free_tag(m, 0, "grad_dye");
  mesh_recover_by_volume(m, "grad_grad_dye");
  mesh_free_tag(m, mesh_dim(m), "grad_grad_dye");
  double weight[1] = {0.075 / 100.0};
  mesh_size_from_hessian(m, "grad_grad_dye", weight, 0.025, 0.1);
  mesh_free_tag(m, 0, "grad_grad_dye");
}

static void warped_adapt(struct mesh** p_m)
{
  for (unsigned i = 0; i < 3; ++i) {
    unsigned done = mesh_warp_to_limit(*p_m, warp_qual_floor);
    set_size_field(*p_m);
    printf("warp to limit, new size field\n");
    write_vtk_step(*p_m);
    mesh_adapt(p_m, size_floor, good_qual_floor, nsliver_layers, max_ops);
    if (done)
      return;
  }
  fprintf(stderr, "warped_adapt still not done after 3 iters\n");
  abort();
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  mesh_derive_model(m, PI / 4);
  mesh_set_rep(m, MESH_FULL);
  mesh_eval_field(m, 0, "adapt_size", 1, size_fun);
  while (refine_by_size(&m, 0));
  start_vtk_steps("warp");
  mesh_eval_field(m, 0, "dye", 1, dye_fun);
  set_size_field(m);
  write_vtk_step(m);
  mesh_adapt(&m, size_floor, good_qual_floor, nsliver_layers, max_ops);
  for (unsigned i = 0; i < 2; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      mesh_eval_field(m, 0, "warp", 3, warp_fun);
      printf("new warp field\n");
      warped_adapt(&m);
      mesh_free_tag(m, 0, "warp");
    }
    the_rotation = -the_rotation;
  }
  free_mesh(m);
}
