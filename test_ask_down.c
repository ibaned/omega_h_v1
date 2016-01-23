#include <stdio.h>
#include "eval_field.h"
#include "mesh.h"
#include "vtk.h"
#include "graph.h"

int main()
{
  printf("Started\n");
  unsigned dim = 3;
  struct mesh* m= new_box_mesh(dim);
  unsigned const* out = mesh_ask_down(m, 3, 0);
  unsigned nverts = mesh_count(m, dim);
  for (unsigned i = 0; i < nverts; ++i) {
    printf(" %u(%u)",i, out[i]);
    printf("\n");
  }
  free_mesh(m);
}
