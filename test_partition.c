#include <assert.h>

#include "comm.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  comm_init();
  assert(argc == 3);
  struct mesh* m = read_and_partition_serial_mesh(argv[1]);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  comm_fini();
}
