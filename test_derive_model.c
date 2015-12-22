#include <assert.h>
#include <stdlib.h>

#include "algebra.h"
#include "ints.h"
#include "loop.h"
#include "derive_model.h"
#include "mesh.h"
#include "subset.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 3);
  unsigned mesh_dim = (unsigned) atoi(argv[1]);
  assert(1 <= mesh_dim);
  assert(mesh_dim <= 3);
  unsigned out_dim = (unsigned) atoi(argv[2]);
  assert(1 <= out_dim);
  assert(out_dim <= 3);
  struct mesh* m = new_box_mesh(mesh_dim);
  mesh_derive_model(m, PI / 4);
  unsigned* offsets = uints_linear(mesh_count(m, out_dim) + 1, 1);
  struct mesh* sm = subset_mesh(m, out_dim, offsets);
  loop_free(offsets);
  free_mesh(m);
  write_vtu(sm, "out.vtu");
  free_mesh(sm);
}
