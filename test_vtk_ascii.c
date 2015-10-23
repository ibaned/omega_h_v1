#include "mesh.h"
#include "vtk.h"
#include <assert.h>

int main(int argc, char** argv)
{
  assert(argc == 3);
  struct mesh* m = read_vtu(argv[1]);
  write_vtu_opts(m, argv[2], VTK_ASCII);
  free_mesh(m);
}
