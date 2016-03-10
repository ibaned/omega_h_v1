#include <assert.h>

#include "arrays.hpp"
#include "include/omega_h.hpp"
#include "mesh.hpp"
#include "coarsen.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  assert(argc == 3);
  osh_init(&argc, &argv);
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_add_tag(m, 0, TAG_F64, "adapt_size", 1,
      doubles_filled(mesh_count(m, 0), 2.0));
  coarsen_by_size(m, 0.1, 0.7);
  mesh_free_tag(m, 0, "adapt_size");
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}

