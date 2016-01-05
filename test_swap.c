#include <assert.h>
#include <stdio.h>

#include "ints.h"
#include "mesh.h"
#include "swap_common.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 3);
  struct mesh* m = read_vtu(argv[1]);
  unsigned* candidates = uints_filled(mesh_count(m, 1), 1);
  unsigned did = swap_common(&m, candidates);
  if (did)
    printf("swapped!\n");
  loop_free(candidates);
  write_vtu(m, argv[2]);
  free_mesh(m);
}
