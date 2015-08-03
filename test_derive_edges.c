#include "derive_edges.h"
#include "tables.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
  unsigned elem_dim = 3;
  struct derived_edges de = derive_edges(
      elem_dim,
      the_box_nelems[elem_dim],
      the_box_nverts[elem_dim],
      the_box_conns[elem_dim]);
  printf("%u edges\n", de.nedge);
  unsigned const* ev = de.edge_verts;
  for (unsigned i = 0; i < de.nedge; ++i) {
    printf("[%u] = %u %u\n", i, ev[0], ev[1]);
    ev += 2;
  }
  free(de.edge_verts);
}
