#include "mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  (void) argc;
  struct mesh* m = read_vtu(argv[1]);
  write_vtu(m, argv[2]);
  free_mesh(m);
}
