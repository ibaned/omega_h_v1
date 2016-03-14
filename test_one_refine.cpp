#include <cassert>

#include "include/omega_h.h"
#include "mesh.hpp"
#include "refine.hpp"
#include "vtk_io.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  assert(argc == 3);
  osh_init(&argc, &argv);
  struct mesh* m = read_mesh_vtk(argv[1]);
  uniformly_refine(m);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}
