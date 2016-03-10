#include <assert.h>

#include "gmsh_io.hpp"
#include "mesh.hpp"
#include "include/omega_h.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  assert(argc == 3);
  osh_init();
  struct mesh* m = read_msh(argv[1]);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}
