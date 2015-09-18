#include "vtk.h"
#include "mesh.h"
#include "mark.h"
#include "subset.h"
#include "loop.h"
#include <assert.h>

int main(int argc, char** argv)
{
  assert(argc == 3);
  struct mesh* m = read_vtk(argv[1]);
  unsigned* pb = mesh_mark_part_boundary(m);
  struct mesh* sm = subset_mesh(m, mesh_dim(m) - 1, pb);
  free_mesh(m);
  loop_free(pb);
  write_vtk(sm, argv[2]);
  free_mesh(sm);
}
