#include <assert.h>

#include "comm.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  comm_init();
  assert(argc == 3);
  struct mesh* m = read_and_partition_serial_mesh(argv[1]);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  comm_fini();
}
