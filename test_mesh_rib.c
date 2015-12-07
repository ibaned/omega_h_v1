#include <assert.h>

#include "bcast.h"
#include "comm.h"
#include "mesh.h"
#include "parallel_inertial_bisect.h"
#include "parallel_mesh.h"
#include "refine_by_size.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  comm_init();
  assert(argc == 3);
  struct mesh* m = 0;
  if (comm_rank() == 0)
    m = read_vtu(argv[1]);
  m = bcast_mesh_metadata(m);
  mesh_number_simply(m, 0);
  balance_mesh_inertial(&m);
  mesh_global_renumber(m, 0);
  write_parallel_vtu(m, argv[2]);
  free_mesh(m);
  comm_fini();
}
