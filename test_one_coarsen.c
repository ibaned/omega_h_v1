#include <assert.h>

#include "arrays.h"
#include "comm.h"
#include "mesh.h"
#include "coarsen.h"
#include "vtk_io.h"

int main(int argc, char** argv)
{
  assert(argc == 3);
  comm_init();
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_add_tag(m, 0, TAG_F64, "adapt_size", 1,
      doubles_filled(mesh_count(m, 0), 2.0));
  coarsen_by_size(m, 0.1, 0.7);
  mesh_free_tag(m, 0, "adapt_size");
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  comm_fini();
}

