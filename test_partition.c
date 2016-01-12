#include <assert.h>

#include "comm.h"
#include "mesh.h"
#include "parallel_inertial_bisect.h"
#include "parallel_mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  comm_init();
  assert(argc == 3);
  struct mesh* m = 0;
  if (comm_rank() == 0) {
    comm_use(comm_self());
    m = read_mesh_vtk(argv[1]);
    assert(!mesh_is_parallel(m));
    mesh_make_parallel(m);
    comm_use(comm_world());
  }
  mesh_partition_out(&m, comm_size());
  balance_mesh_inertial(&m);
  for (unsigned d = 0; d <= mesh_dim(m); ++d)
    if (mesh_has_dim(m, d))
      mesh_global_renumber(m, d);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  comm_fini();
}

