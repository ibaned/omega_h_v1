#include "include/omega_h.h"
#include "comm.hpp"
#include "vtk_io.hpp"
#include "parallel_mesh.hpp"
#include "parallel_inertial_bisect.hpp"
#include "mesh.hpp"

#include <cstdio>
#include <cassert>

using namespace omega_h;

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  unsigned world_size = comm_size();
  unsigned subcomm_size = 0;
  struct mesh* m = 0;
  for (subcomm_size = 1; subcomm_size <= world_size; subcomm_size *= 2) {
    unsigned group = comm_rank() / subcomm_size;
    struct comm* subcomm = comm_split(comm_using(), group, comm_rank() % subcomm_size);
    comm_use(subcomm);
    if (!group) {
      if (comm_size() == 1) {
        m = read_mesh_vtk(argv[1]);
        mesh_make_parallel(m);
      } else {
        mesh_partition_out(&m, m != 0);
        balance_mesh_inertial(m);
      }
      char fname[256];
      sprintf(fname, "out_%u.pvtu", subcomm_size);
      write_mesh_vtk(m, fname);
    }
    comm_use(comm_world());
    comm_free(subcomm);
  }
  assert(m);
  free_mesh(m);
  osh_fini();
}
