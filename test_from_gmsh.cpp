#include <cassert>

#include "gmsh_io.hpp"
#include "mesh.hpp"
#include "include/omega_h.h"
#include "vtk_io.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 3);
  struct mesh* m = read_msh(argv[1]);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}
