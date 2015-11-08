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
#ifdef __CUDACC__
  unsigned const * in = (unsigned const *)
		  loop_cuda_to_device( the_box_conns[dim] , sizeof(unsigned) * nelems * nverts);
#else
  unsigned const *in = the_box_conns[dim];
#endif
  up_from_down(
      dim,
      0,
      nelems,
      nverts,
      in,
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
  loop_free(offsets);
  loop_free(adj);
  loop_free(directions);
}
