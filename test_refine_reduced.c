#include "refine_reduced.h"
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
  struct rv_mesh m = new_box_rv_mesh(3);
  char fname[64];
  for (unsigned it = 0; 1; ++it) {
    struct rv_mesh out = refine_reduced(m, simple);
    printf("%u elements, %u vertices\n", out.nelems, out.nverts);
    double* quals = element_qualities(
        out.elem_dim,
        out.nelems,
        out.verts_of_elems,
        out.coords);
    double minqual = doubles_min(quals, m.nelems);
    free(quals);
    printf("min quality %f\n", minqual);
    if (out.nelems == m.nelems) {
      free_rv_mesh(out);
      break;
    }
    free_rv_mesh(m);
    m = out;
    sprintf(fname, "out_%u.vtu", it);
    write_vtk(fname, m.elem_dim, m.nelems, m.nverts,
        m.verts_of_elems, m.coords);
  }
  free_rv_mesh(m);
}
