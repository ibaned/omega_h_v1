#include <cmath>
#include <cstdio>
#include <time.h>

#include "algebra.hpp"
#include "derive_model.hpp"
#include "eval_field.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "include/omega_h.h"
#include "refine.hpp"
#include "vtk_io.hpp"

#ifdef LOOP_CUDA_HPP
static double get_time(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  double s = ts.tv_sec;
  double ns = ts.tv_nsec;
  s += (ns / 1e9);
  return s;
}

static void trigger_cuda_init(void)
{
  double t0 = get_time();
  int* p;
  cudaMalloc(&p, sizeof(int));
  cudaFree(p);
  double t1 = get_time();
  printf("CUDA takes %f seconds to initialize\n", t1 - t0);
}
#endif

static void sinusoid(double const* x, double* size)
{
  double s = cos(x[0] * 8.0 * PI) / 4.0 + 1.0 / 2.0;
  double d = fabs(x[1] - s);
  double fudge = 1.4;
  *size = sqrt(2 * (1e-7 + d * 1e-5)) * fudge;
}

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
#ifdef LOOP_CUDA_HPP
  trigger_cuda_init();
#endif
  struct mesh* m = new_box_mesh(2);
  mesh_derive_model(m, PI / 4);
  mesh_set_rep(m, MESH_FULL);
  unsigned did_refine = 0;
#ifdef LOOP_CUDA_HPP
  double t0 = get_time();
#endif
  do {
    mesh_eval_field(m, 0, "adapt_size", 1, sinusoid);
    did_refine = refine_by_size(m, 0.0);
  //printf("%u triangles\n", mesh_count(m, 2));
    mesh_free_tag(m, 0, "adapt_size");
  } while (did_refine);
#ifdef LOOP_CUDA_HPP
  double t1 = get_time();
  printf("refinement time %f seconds, %u final triangles\n",
      t1 - t0, mesh_count(m, 2));
#endif
  free_mesh(m);
  osh_fini();
}
