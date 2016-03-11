#include <cassert>

#include "include/omega_h.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 3);
  struct mesh* m = read_and_partition_serial_mesh(argv[1]);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}
