#include <cassert>
#include <cstdlib>

#include "ghost_mesh.hpp"
#include "int_casts.hpp"
#include "mesh.hpp"
#include "include/omega_h.h"
#include "vtk_io.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 4);
  struct mesh* m = 0;
  m = read_mesh_vtk(argv[1]);
  unsigned nlayers = U(atoi(argv[2]));
  assert(nlayers <= 10);
  mesh_ensure_ghosting(m, nlayers);
  write_mesh_vtk(m, argv[3]);
  free_mesh(m);
  osh_fini();
}
