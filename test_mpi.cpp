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
  if (!comm_rank())
    printf("#elements: total %lu, max %u\n",
        total, max);
}

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 3);
  const char* filename = argv[1];
  unsigned const step_factor = static_cast<unsigned>(atoi(argv[2]));
  unsigned world_size = comm_size();
  unsigned subcomm_size = 0;
  struct mesh* m = 0;
  double total_time = 0;
  double total_balance_time = 0;
  double total_refine_time = 0;
  unsigned long target_nelems = 0;
  unsigned starting_nelems = 0;
  if (comm_rank() == 0) {
    comm_use(comm_self());
    double q = get_time();
    m = read_mesh_vtk(filename);
    mesh_make_parallel(m);
    double r = get_time();
    printf("setup time %.4e seconds\n", r - q);
    print_stats(m);
    comm_use(comm_world());
    starting_nelems = mesh_count(m, mesh_dim(m));
  }
  target_nelems = comm_bcast_uint(starting_nelems);
  for (subcomm_size = step_factor;
       subcomm_size <= world_size;
       subcomm_size *= step_factor) {
    double a = get_time();
    target_nelems *= step_factor;
    if (!comm_rank())
      printf("subcomm size %u, target %lu\n", subcomm_size, target_nelems);
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
        printf("balance time %.4e seconds\n", balance_time);
      print_stats(m);
      do {
        double e = get_time();
        uniformly_refine(m, 0.2);
        double f = get_time();
        double refine_time = comm_max_double(f - e);
        double minqual = comm_min_double(mesh_min_quality(m));
        total_refine_time += refine_time;
        if (!comm_rank())
          printf("refine time %.4e seconds, minqual %e\n", refine_time, minqual);
        print_stats(m);
      } while (comm_add_ulong(mesh_count(m, mesh_dim(m))) < target_nelems);
    }
    comm_use(comm_world());
    comm_free(subcomm);
    double b = get_time();
    double step_time = comm_max_double(b - a);
    total_time += step_time;
    if (!comm_rank())
      printf("step time %.4e seconds\n", step_time);
  }
  if (!comm_rank()) {
    printf("total time %.4e seconds\n", total_time);
    printf("total refine time %.4e seconds\n", total_refine_time);
    printf("total balance time %.4e seconds\n", total_balance_time);
  }
  assert(m);
  free_mesh(m);
  osh_fini();
}
