#include "mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  (void) argc;
  struct mesh* m = read_mesh_vtk(argv[1]);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
}
