#include <assert.h>
#include <stdlib.h>

#include "ghost_mesh.h"
#include "mesh.h"
#include "include/omega_h.h"
#include "vtk_io.h"

int main(int argc, char** argv)
{
  osh_init();
  assert(argc == 4);
  struct mesh* m = 0;
  m = read_mesh_vtk(argv[1]);
  unsigned nlayers = (unsigned) atoi(argv[2]);
  assert(nlayers <= 10);
  mesh_ensure_ghosting(m, nlayers);
  write_mesh_vtk(m, argv[3]);
  free_mesh(m);
  osh_fini();
}
