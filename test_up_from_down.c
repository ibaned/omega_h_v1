#include <stdio.h>

#include "loop.h"
#include "tables.h"
#include "up_from_down.h"

int main()
{
  unsigned dim = 3;
  unsigned nelems = the_box_nelems[dim];
  unsigned nverts = the_box_nverts[dim];
  unsigned* offsets;
  unsigned* adj;
  unsigned* directions;
  up_from_down(
      dim,
      0,
      nelems,
      nverts,
      the_box_conns[dim],
      &offsets,
      &adj,
      &directions);
  for (unsigned i = 0; i < nverts; ++i) {
    printf("[%u] =", i);
    unsigned first_adj = offsets[i];
    unsigned last_adj = offsets[i + 1];
    for (unsigned j = first_adj; j < last_adj; ++j)
      printf(" %u(%u)", adj[j], directions[j]);
    printf("\n");
  }
  LOOP_FREE(offsets);
  LOOP_FREE(adj);
  LOOP_FREE(directions);
}
