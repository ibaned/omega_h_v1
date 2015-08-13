#include "refine_by_size.h"
#include "tables.h"
#include "algebra.h"
#include "vtk.h"
#include "element_qualities.h"
#include "doubles.h"
#include <stdlib.h>
#include <stdio.h>

static double simple(double const x[])
{
  double coarse = 0.5;
  double fine = 0.025;
  double radius = vector_norm(x, 3);
  double d = fabs(radius - 0.5);
  return coarse * d + fine * (1 - d);
}

int main()
{
  unsigned const elem_dim = 3;
  unsigned nelems;
  unsigned nverts;
  unsigned* verts_of_elems;
  double* coords;
  get_box_copy(elem_dim, &nelems, &nverts, &verts_of_elems, &coords);
  char fname[64];
  for (unsigned it = 0; 1; ++it) {
    if (!refine_by_size(elem_dim,
          &nelems, &nverts, &verts_of_elems, &coords,
          simple))
      break;
    printf("%u elements, %u vertices\n", nelems, nverts);
    double* quals = element_qualities(elem_dim, nelems, verts_of_elems, coords);
    double minqual = doubles_min(quals, nelems);
    free(quals);
    printf("min quality %f\n", minqual);
    sprintf(fname, "out_%u.vtu", it);
    write_vtk(fname, elem_dim, nelems, nverts, verts_of_elems, coords, 0);
  }
  free(verts_of_elems);
  free(coords);
}
