#include "refine_topology.h"
#include "tables.h"
#include <stdio.h>
#include <stdlib.h>

int main()
{
  unsigned elem_dim = 3;
  unsigned split_dim = 3;
  unsigned ent_dim = 3;
  unsigned nelem = 1;
  unsigned const elem_verts[4] = {0,1,2,3};
  unsigned const elem_split_offset[2] = {0,1};
  unsigned const elem_split_vert[1] = {4};
  unsigned const elem_split_direction[1] = {0};
  unsigned nprods;
  unsigned* verts_of_prods;
  refine_topology(
      elem_dim,
      split_dim,
      ent_dim,
      nelem,
      elem_verts,
      elem_split_offset,
      elem_split_vert,
      elem_split_direction,
      &nprods,
      &verts_of_prods);
  unsigned verts_per_prod = the_down_degrees[ent_dim][0];
  printf("%u entities:\n", nprods);
  for (unsigned i = 0; i < nprods; ++i) {
    printf("[%u] =", i);
    for (unsigned j = 0; j < verts_per_prod; ++j)
      printf(" %u", verts_of_prods[i * verts_per_prod + j]);
    printf("\n");
  }
  free(verts_of_prods);
}
