#include "derive_edges.h"
#include "tables.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
  unsigned elem_dim = 3;
  unsigned nedges;
  unsigned* verts_of_edges;
  derive_edges(
      elem_dim,
      the_box_nelems[elem_dim],
      the_box_nverts[elem_dim],
      the_box_conns[elem_dim],
      &nedges,
      &verts_of_edges);
  printf("%u edges\n", nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    unsigned const* ev = verts_of_edges + i * 2;
    printf("[%u] = %u %u\n", i, ev[0], ev[1]);
  }
  free(verts_of_edges);
}
