#include "up_from_down.h"
#include "tables.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
  unsigned dim = 3;
  unsigned nup = the_box_nelems[dim];
  unsigned ndown = the_box_nverts[dim];
  struct up_adj out = up_from_down(
      dim,
      0,
      nup,
      ndown,
      the_box_conns[dim]);
  for (unsigned i = 0; i < ndown; ++i) {
    printf("[%u] =", i);
    unsigned first_up = out.offsets[i];
    unsigned last_up = out.offsets[i + 1];
    for (unsigned j = first_up; j < last_up; ++j)
      printf(" %u(%u)", out.edges[j], out.directions[j]);
    printf("\n");
  }
  free(out.offsets);
  free(out.edges);
  free(out.directions);
}
