#include <assert.h>
#include <stdlib.h>

#include "algebra.h"
#include "derive_model.h"
#include "mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 4);
  struct mesh* m = read_vtu(argv[1]);
  mesh_derive_model(m, PI / atof(argv[2]));
  write_vtu(m, argv[3]);
  free_mesh(m);
}
