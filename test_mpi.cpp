#include "include/omega_h.h"
#include "comm.hpp"
#include "parallel_mesh.hpp"
#include "parallel_inertial_bisect.hpp"
#include "mesh.hpp"
#include "refine.hpp"
#include "derive_model.hpp"
#include "algebra.hpp"

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
  assert(argc == 2);
  unsigned const initial_refines = static_cast<unsigned>(atoi(argv[1]));
  unsigned world_size = comm_size();
  unsigned subcomm_size = 0;
  struct mesh* m = 0;
  double total_time = 0;
  double total_balance_time = 0;
  double total_refine_time = 0;
  for (subcomm_size = 1; subcomm_size <= world_size; subcomm_size *= 2) {
    double a = get_time();
    if (!comm_rank())
      printf("subcomm size %u\n", subcomm_size);
    unsigned group = comm_rank() / subcomm_size;
    struct comm* subcomm = comm_split(comm_using(), group, comm_rank() % subcomm_size);
    comm_use(subcomm);
    if (!group) {
      if (comm_size() == 1) {
        double q = get_time();
        m = new_box_mesh(3);
        mesh_derive_model(m, PI / 4);
        mesh_set_rep(m, MESH_FULL);
        for (unsigned i = 0; i < initial_refines; ++i)
          uniformly_refine(m);
        mesh_make_parallel(m);
        double r = get_time();
        printf("setup time %.4e seconds\n", r - q);
        print_stats(m);
      } else {
        mesh_partition_out(&m, m != 0);
        double c = get_time();
        balance_mesh_inertial(m);
        double d = get_time();
        double balance_time = comm_max_double(d - c);
        total_balance_time += balance_time;
        if (!comm_rank())
          printf("balance time %.4e seconds\n", balance_time);
        print_stats(m);
      }
      double e = get_time();
      uniformly_refine(m);
      double f = get_time();
      double refine_time = comm_max_double(f - e);
      total_refine_time += refine_time;
      if (!comm_rank())
        printf("refine time %.4e seconds\n", refine_time);
      print_stats(m);
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
