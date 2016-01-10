#include <assert.h>
#include <stdlib.h>

#include "algebra.h"
#include "derive_model.h"
#include "mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 4);
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_derive_model(m, PI / atof(argv[2]));
  write_mesh_vtk(m, argv[3]);
  free_mesh(m);
}
