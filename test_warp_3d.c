#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "adapt.h"
#include "algebra.h"
#include "comm.h"
#include "derive_model.h"
#include "doubles.h"
#include "element_field.h"
#include "eval_field.h"
#include "mesh.h"
#include "refine.h"
#include "vtk.h"
#include "warp_to_limit.h"

struct mesh;

static double const warp_qual_floor = 0.2;
static double const good_qual_floor = 0.3;
static double const size_floor = 1. / 3.;
static unsigned const nsliver_layers = 4;
static unsigned const max_ops = 50;

static void size_fun(double const* x, double* s)
{
  (void)x;
  s[0] = 0.2;
}

static double the_rotation = PI / 4;

static void warp_fun(double const* coords, double* v)
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

static void dye_fun(double const* coords, double* v)
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

static void warped_adapt(struct mesh** p_m)
{
  static unsigned const n = 6;
  for (unsigned i = 0; i < n; ++i) {
    printf("\n WARP TO LIMIT %u\n", i);
    unsigned done = mesh_warp_to_limit(*p_m, warp_qual_floor);
    mesh_adapt(*p_m, size_floor, good_qual_floor, nsliver_layers, max_ops);
    write_vtk_step(*p_m);
    if (done)
      return;
  }
  fprintf(stderr, "warped_adapt still not done after %u iters\n", n);
  abort();
}

int main()
{
  comm_init();
  struct mesh* m = new_box_mesh(3);
  mesh_derive_model(m, PI / 4);
  mesh_set_rep(m, MESH_FULL);
  mesh_eval_field(m, 0, "adapt_size", 1, size_fun);
  while (refine_by_size(m, 0));
  start_vtk_steps("warp");
  mesh_eval_field(m, 0, "dye", 1, dye_fun);
  write_vtk_step(m);
  for (unsigned i = 0; i < 2; ++i) {
    printf("\nOUTER DIRECTION %u\n", i);
    for (unsigned j = 0; j < 4; ++j) {
      printf("\nWARP FIELD %u\n", j);
      mesh_eval_field(m, 0, "warp", 3, warp_fun);
      printf("new warp field\n");
      warped_adapt(&m);
      mesh_free_tag(m, 0, "warp");
    }
    the_rotation = -the_rotation;
  }
  free_mesh(m);
  comm_fini();
}
