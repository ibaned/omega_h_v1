#include <assert.h>
#include <stdlib.h>

#include "bcast.h"
#include "comm.h"
#include "ghost_mesh.h"
#include "mesh.h"
#include "parallel_inertial_bisect.h"
#include "parallel_mesh.h"
#include "refine_by_size.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  comm_init();
  assert(argc == 4);
  struct mesh* m = 0;
  if (comm_rank() == 0) {
    m = read_vtu(argv[1]);
    mesh_make_parallel(m);
  }
  m = bcast_mesh_metadata(m);
  balance_mesh_inertial(&m);
  mesh_global_renumber(m, 0);
  unsigned nlayers = (unsigned) atoi(argv[2]);
  assert(nlayers <= 10);
  ghost_mesh(&m, nlayers);
  write_parallel_vtu(m, argv[3]);
  free_mesh(m);
  comm_fini();
}
