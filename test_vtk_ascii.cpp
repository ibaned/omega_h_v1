#include <assert.h>

#include "comm.h"
#include "mesh.h"
#include "vtk_io.h"

int main(int argc, char** argv)
{
  comm_init();
  assert(argc == 3);
  struct mesh* m = read_mesh_vtk(argv[1]);
  write_mesh_vtk_opts(m, argv[2], VTK_ASCII);
  free_mesh(m);
  comm_fini();
}
