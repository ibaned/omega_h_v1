#include <cassert>

#include "arrays.hpp"
#include "include/omega_h.h"
#include "mesh.hpp"
#include "coarsen.hpp"
#include "vtk_io.hpp"
#include "adapt.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  assert(argc == 2);
  osh_init(&argc, &argv);
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_add_tag(m, 0, TAG_F64, "adapt_size", 1,
      filled_array(mesh_count(m, 0), 0.1));
  start_vtk_steps("imr");
  for (unsigned i = 0; i < 1; ++i) {
    mesh_adapt(m, 0.3, 0.3, 4, 50);
    write_vtk_step(m);
  }
  free_mesh(m);
  osh_fini();
}


