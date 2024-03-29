#include <cassert>
#include <cstdio>

#include "bcast.hpp"
#include "comm.hpp"
#include "include/omega_h.h"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "migrate_mesh.hpp"
#include "parallel_mesh.hpp"
#include "subset.hpp"
#include "vtk_io.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  char const* path;
  if (argc == 2)
    path = argv[1];
  else
    path = ".";
  assert(comm_size() == 2);
  struct mesh* m = 0;
  if (comm_rank() == 0) {
    m = new_box_mesh(2);
    mesh_count(m, 1); //trigger edge creation
  }
  m = bcast_mesh_metadata(m, comm_rank() == 0);
  mesh_make_parallel(m);
  char file[128];
  sprintf(file, "%s/before.pvtu", path);
  write_mesh_vtk(m, file);
  unsigned n = 1;
  if (comm_rank() == 0) {
    unsigned recvd_elem_ranks[1] = {0};
    unsigned recvd_elem_ids[1] = {0};
    migrate_mesh(m, n, recvd_elem_ranks, recvd_elem_ids);
  } else {
    unsigned recvd_elem_ranks[1] = {0};
    unsigned recvd_elem_ids[1] = {1};
    migrate_mesh(m, n, recvd_elem_ranks, recvd_elem_ids);
  }
  sprintf(file, "%s/after.pvtu", path);
  write_mesh_vtk(m, file);
  free_mesh(m);
  osh_fini();
}
