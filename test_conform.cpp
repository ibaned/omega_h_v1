#include <cassert>
#include <cstdio>

#include "arrays.hpp"
#include "bcast.hpp"
#include "comm.hpp"
#include "global.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "migrate_mesh.hpp"
#include "include/omega_h.hpp"
#include "parallel_mesh.hpp"
#include "vtk_io.hpp"

static struct mesh* make_2_tri_parallel(void)
{
  struct mesh* m = 0;
  assert(comm_size() == 2);
  if (comm_rank() == 0) {
    m = new_box_mesh(2);
  }
  m = bcast_mesh_metadata(m);
  mesh_make_parallel(m);
  if (comm_rank() == 0) {
    unsigned n = 1;
    unsigned recvd_elem_ranks[1] = {0};
    unsigned recvd_elem_ids[1] = {0};
    migrate_mesh(m, n, recvd_elem_ranks, recvd_elem_ids);
  } else {
    unsigned n = 1;
    unsigned recvd_elem_ranks[1] = {0};
    unsigned recvd_elem_ids[1] = {1};
    migrate_mesh(m, n, recvd_elem_ranks, recvd_elem_ids);
  }
  return m;
}

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  char const* path;
  if (argc == 2)
    path = argv[1];
  else
    path = ".";
  struct mesh* m = make_2_tri_parallel();
  double* data = doubles_filled(3, static_cast<double>(comm_rank() + 1));
  mesh_add_tag(m, 0, TAG_F64, "field", 1, data);
  char file[128];
  sprintf(file, "%s/one.pvtu", path);
  write_mesh_vtk(m, file);
  mesh_accumulate_tag(m, 0, "field");
  sprintf(file, "%s/two.pvtu", path);
  write_mesh_vtk(m, file);
  mesh_conform_tag(m, 0, "field");
  sprintf(file, "%s/three.pvtu", path);
  write_mesh_vtk(m, file);
  free_mesh(m);
  osh_fini();
}
