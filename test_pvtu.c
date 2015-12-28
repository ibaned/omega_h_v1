#include <assert.h>

#include "mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 3);
  struct mesh* m = read_vtu(argv[1]);
  write_pvtu(m, argv[2], 100);
  free_mesh(m);
}
