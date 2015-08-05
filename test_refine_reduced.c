#include "refine_reduced.h"
#include "tables.h"
#include <stdlib.h>

static double sf(double const x[])
{
  (void) x;
  return 1.0;
}

int main()
{
  unsigned dim = 2;
  unsigned nelems;
  unsigned nverts;
  unsigned* verts_of_elems;
  double* coords;
  refine_reduced(
      dim,
      the_box_nelems[dim],
      the_box_nverts[dim],
      the_box_conns[dim],
      the_box_coords[dim],
      sf,
      &nelems,
      &nverts,
      &verts_of_elems,
      &coords);
  free(verts_of_elems);
  free(coords);
}
