#include <assert.h>
#include <stdlib.h>

#include "mesh.h"
#include "include/omega_h.h"
#include "size.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  osh_init();
  assert(argc == 3);
  struct mesh* m = 0;
  m = read_mesh_vtk(argv[1]);
  mesh_identity_size_field(m, "identity_size");
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}

