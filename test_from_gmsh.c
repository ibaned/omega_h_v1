#include <assert.h>

#include "gmsh_io.h"
#include "mesh.h"
#include "refine_by_class.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 3);
  struct mesh* m = read_msh(argv[1]);
//refine_by_class(&m);
  write_vtu(m, argv[2]);
  free_mesh(m);
}
