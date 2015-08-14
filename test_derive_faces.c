#include "tables.h"
#include "up_from_down.h"
#include "reflect_down.h"
#include "bridge_graph.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
  unsigned elem_dim = 3;
  unsigned nelems;
  unsigned nverts;
  unsigned* verts_of_elems;
  get_box_copy(elem_dim, &nelems, &nverts, &verts_of_elems, 0);
  unsigned* elems_of_verts_offsets;
  unsigned* elems_of_verts;
  up_from_down(elem_dim, 0, nelems, nverts, verts_of_elems,
      &elems_of_verts_offsets, &elems_of_verts, 0);
  unsigned* elems_of_elems = get_dual(elem_dim, nelems,
      verts_of_elems, elems_of_verts_offsets, elems_of_verts);
  free(verts_of_elems);
  free(elems_of_verts_offsets);
  free(elems_of_verts);
  for (unsigned i = 0; i < nelems; ++i) {
    printf("%u:", i);
    for (unsigned j = 0; j < 4; ++j) {
      unsigned ei = elems_of_elems[i * 4 + j];
      if (ei == INVALID)
        printf(" X");
      else
        printf(" %u", ei);
    }
    printf("\n");
  }
  free(elems_of_elems);
}
