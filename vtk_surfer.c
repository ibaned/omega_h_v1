#include "vtk.h"
#include "mesh.h"
#include "mark.h"
#include "subset.h"
#include "loop.h"
#include "ints.h"
#include <assert.h>

int main(int argc, char** argv)
{
  assert(argc == 3);
  struct mesh* m = read_vtk(argv[1]);
  unsigned d = mesh_dim(m);
  unsigned ns = mesh_count(m, d - 1);
  unsigned* pb = mesh_mark_part_boundary(m);
  unsigned* off = ints_exscan(pb, ns);
  loop_free(pb);
  struct mesh* sm = subset_mesh(m, d - 1, off);
  free_mesh(m);
  loop_free(off);
  write_vtk(sm, argv[2]);
  free_mesh(sm);
}
