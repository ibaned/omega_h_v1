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
#ifdef __CUDACC__
  unsigned t_size  = the_down_degrees[dim][0] * nelems;
  unsigned*t_off = offsets;
  offsets = (unsigned*)
    loop_cuda_to_host( offsets , sizeof(unsigned)* nverts);
  unsigned*t_adj = adj;
  adj = (unsigned*)
    loop_cuda_to_host( adj , sizeof(unsigned)* t_size);
  unsigned*t_directions = directions;
  directions = (unsigned*)
    loop_cuda_to_host( directions , sizeof(unsigned)* t_size);
#endif
  for (unsigned i = 0; i < nverts; ++i) {
    printf("[%u] =", i);
    unsigned first_adj = offsets[i];
    unsigned last_adj = offsets[i + 1];
    for (unsigned j = first_adj; j < last_adj; ++j)
      printf(" %u(%u)", adj[j], directions[j]);
    printf("\n");
  }
#ifdef __CUDACC__
  loop_free(t_off);
  loop_free(t_adj);
  loop_free(t_directions);
  free(offsets);free(adj);free(directions);
#else
  loop_free(offsets);
  loop_free(adj);
  loop_free(directions);
#endif
}
