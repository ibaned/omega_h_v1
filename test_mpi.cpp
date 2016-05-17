#include "include/omega_h.h"
#include "comm.hpp"
#include "parallel_mesh.hpp"
#include "parallel_inertial_bisect.hpp"
#include "mesh.hpp"
#include "refine.hpp"
#include "derive_model.hpp"
#include "algebra.hpp"
#include "vtk_io.hpp"
#include "quality.hpp"

#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>

using namespace omega_h;

static double get_time(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double t = static_cast<double>(tv.tv_usec);
  t /= 1e6;
  t += static_cast<double>(tv.tv_sec);
  return t;
}

static void print_stats(struct mesh* m)
{
  unsigned nelems = mesh_count(m, mesh_dim(m));
  unsigned long total = comm_add_ulong(nelems);
  unsigned max = comm_max_uint(nelems);
  double minqual = comm_min_double(mesh_min_quality(m));
  if (!comm_rank())
    printf("#elements: total %lu, max %u | min quality %e\n",
        total, max, minqual);
}

LOOP_KERNEL(decrease_vert_size_field,
    double* sf)
  sf[i] /= cbrt(2);
}

static void decrease_size_field(struct mesh* m)
{
  if (!comm_rank())
    printf("decreasing size field by 2^(1/3)\n");
  double* sf = mesh_find_tag(m, 0, "adapt_size")->d.f64;
  unsigned nverts = mesh_count(m, 0);
  LOOP_EXEC(decrease_vert_size_field, nverts, sf);
}

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 6);
  const char* filename = argv[1];
  unsigned const initial_refs = static_cast<unsigned>(atoi(argv[2]));
  unsigned const step_factor = static_cast<unsigned>(atoi(argv[3]));
  double const qual_floor = atof(argv[4]);
  const char* outfile = argv[5];
  unsigned world_size = comm_size();
  unsigned subcomm_size = 0;
  struct mesh* m = 0;
  double total_time = 0;
  double total_balance_time = 0;
  double total_refine_time = 0;
  if (comm_rank() == 0) {
    comm_use(comm_self());
    double q = get_time();
    m = read_mesh_vtk(filename);
    print_stats(m);
    mesh_identity_size_field(m, "adapt_size");
    for (unsigned i = 0; i < initial_refs; ++i) {
      decrease_size_field(m);
      while (refine_by_size(m, qual_floor)) {
        print_stats(m);
      }
    }
    mesh_make_parallel(m);
    double r = get_time();
    printf("setup time %f seconds\n", r - q);
    comm_use(comm_world());
  }
  for (subcomm_size = step_factor;
       subcomm_size <= world_size;
       subcomm_size *= step_factor) {
    double a = get_time();
    if (!comm_rank())
      printf("subcomm size %u\n", subcomm_size);
    unsigned group = comm_rank() / subcomm_size;
    struct comm* subcomm = comm_split(comm_using(), group, comm_rank() % subcomm_size);
    comm_use(subcomm);
    if (!group) {
      mesh_partition_out(&m, m != 0);
      double c = get_time();
      balance_mesh_inertial(m);
      double d = get_time();
      double balance_time = comm_max_double(d - c);
      total_balance_time += balance_time;
      if (!comm_rank())
        printf("balance time %f seconds\n", balance_time);
      print_stats(m);
      decrease_size_field(m);
      while (1) {
        double e = get_time();
        unsigned did = refine_by_size(m, qual_floor);
        double f = get_time();
        print_stats(m);
        double refine_time = comm_max_double(f - e);
        total_refine_time += refine_time;
        if (!comm_rank())
          printf("refine time %f seconds\n", refine_time);
        if (!did)
          break;
      }
    }
    comm_use(comm_world());
    comm_free(subcomm);
    double b = get_time();
    double step_time = comm_max_double(b - a);
    total_time += step_time;
    if (!comm_rank())
      printf("step time %f seconds\n", step_time);
  }
  if (!comm_rank()) {
    printf("total time %f seconds\n", total_time);
    printf("total refine time %f seconds\n", total_refine_time);
    printf("total balance time %f seconds\n", total_balance_time);
  }
  write_mesh_vtk(m, outfile);
  assert(m);
  free_mesh(m);
  osh_fini();
}
