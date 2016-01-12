#include <stdio.h>
#include "eval_field.h"
#include "mesh.h"
#include "refine_by_size.h"
#include "vtk.h"
#include "graph.h"

int main()
{
  printf("Started\n");
  unsigned dim = 3;
  struct mesh* m= new_box_mesh(dim);
  struct const_up* a = mesh_ask_up(m, 0, 3);
  unsigned nverts = mesh_count(m, dim);
  for (unsigned i = 0; i < nverts; ++i) {
    printf("[%u] =", i);
    unsigned first_adj = a->offsets[i];
    unsigned last_adj = a->offsets[i + 1];
    for (unsigned j = first_adj; j < last_adj; ++j)
      printf(" %u(%u,%u)",j, a->adj[j] , a->directions[j]);
    printf("\n");
  }
  free_mesh(m);
}
