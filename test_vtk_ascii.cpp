#include <cassert>

#include "include/omega_h.h"
#include "mesh.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 3);
  struct mesh* m = read_mesh_vtk(argv[1]);
  write_mesh_vtk_opts(m, argv[2], VTK_ASCII);
  free_mesh(m);
  osh_fini();
}
