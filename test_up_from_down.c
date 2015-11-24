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
  unsigned const* in = (unsigned const*)
      loop_to_device(the_box_conns[dim], sizeof(unsigned) * nelems * nverts);
  up_from_down(
      dim,
      0,
      nelems,
      nverts,
      in,
      &offsets,
      &adj,
      &directions);
  unsigned t_size  = the_down_degrees[dim][0] * nelems;
  unsigned* t_off = offsets;
  offsets = (unsigned*)
    loop_to_host(offsets, (sizeof(unsigned) * nverts) + 1);
  unsigned* t_adj = adj;
  adj = (unsigned*)
    loop_to_host(adj, sizeof(unsigned) * t_size);
  unsigned* t_directions = directions;
  directions = (unsigned*)
    loop_to_host(directions, sizeof(unsigned) * t_size);
  for (unsigned i = 0; i < nverts; ++i) {
    printf("[%u] =", i);
    unsigned first_adj = offsets[i];
    unsigned last_adj = offsets[i + 1];
    for (unsigned j = first_adj; j < last_adj; ++j)
      printf(" %u(%u)", adj[j], directions[j]);
    printf("\n");
  }
  loop_free(t_off);
  loop_free(t_adj);
  loop_free(t_directions);
  loop_host_free(offsets);
  loop_host_free(adj);
  loop_host_free(directions);
}
