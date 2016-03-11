#include <cassert>

#include "include/omega_h.hpp"
#include "mesh.hpp"
#include "swap.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 3);
  struct mesh* m = read_mesh_vtk(argv[1]);
  swap_slivers(m, 0.3, 4);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}


