#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "adapt.hpp"
#include "algebra.hpp"
#include "comm.hpp"
#include "derive_model.hpp"
#include "doubles.hpp"
#include "element_field.hpp"
#include "eval_field.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "parallel_inertial_bisect.hpp"
#include "refine.hpp"
#include "vtk_io.hpp"
#include "warp_to_limit.hpp"

using namespace omega_h;

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

static unsigned global_nelems_for_mass;

static void mass_fun(double const* coords, double* v)
{
  if (coords[0] > 0.5)
    *v = 1.0 / global_nelems_for_mass;
  else
    *v = 0;
}

static void pointwise_fun(double const* coords, double* v)
{
  copy_vector(coords, v, 3);
}

static void warped_adapt(struct mesh* m)
{
  static unsigned const n = 6;
  for (unsigned i = 0; i < n; ++i) {
    if (!comm_rank())
      printf("\n WARP TO LIMIT %u\n", i);
    unsigned done = mesh_warp_to_limit(m, warp_qual_floor);
    mesh_adapt(m, size_floor, good_qual_floor, nsliver_layers, max_ops);
    write_vtk_step(m);
    if (done)
      return;
  }
  fprintf(stderr, "warped_adapt still not done after %u iters\n", n);
  abort();
}

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  char const* path;
  if (argc == 2)
    path = argv[1];
  else
    path = ".";
  struct mesh* m = 0;
  if (comm_rank() == 0) {
    int was_parallel = (comm_size() > 1);
    comm_use(comm_self());
    m = new_box_mesh(3);
    mesh_derive_model(m, PI / 4);
    mesh_set_rep(m, MESH_FULL);
    mesh_eval_field(m, 0, "adapt_size", 1, size_fun);
    while (refine_by_size(m, 0));
    if (was_parallel)
      mesh_make_parallel(m);
    comm_use(comm_world());
  }
  if (comm_size() > 1) {
    mesh_partition_out(&m, comm_rank() == 0);
    balance_mesh_inertial(m);
  }
  char prefix[128];
  sprintf(prefix, "%s/warp", path);
  start_vtk_steps(prefix);
  mesh_eval_field(m, 0, "dye", 1, dye_fun);
  {
    mesh_interp_to_elems(m, "coordinates");
    global_nelems_for_mass = mesh_count(m, mesh_dim(m));
    //set mass field to test conservative transfer
    mesh_eval_field2(m, mesh_dim(m), "mass", 1,
        OSH_TRANSFER_CONSERVE, mass_fun);
    //also now test pointwise transfer
    mesh_eval_field2(m, mesh_dim(m), "pointwise", 3,
        OSH_TRANSFER_POINTWISE, pointwise_fun);
    mesh_free_tag(m, mesh_dim(m), "coordinates");
  }
  write_vtk_step(m);
  for (unsigned i = 0; i < 2; ++i) {
    if (!comm_rank())
      printf("\nOUTER DIRECTION %u\n", i);
    for (unsigned j = 0; j < 4; ++j) {
      if (!comm_rank())
        printf("\nWARP FIELD %u\n", j);
      mesh_eval_field(m, 0, "warp", 3, warp_fun);
      if (!comm_rank())
        printf("new warp field\n");
      warped_adapt(m);
      mesh_free_tag(m, 0, "warp");
    }
    the_rotation = -the_rotation;
  }
  free_mesh(m);
  osh_fini();
}
