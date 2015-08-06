#include "refine_reduced.h"
#include "tables.h"
#include "vtk.h"
#include <stdlib.h>
#include <stdio.h>

static double simple(double const x[])
{
  return 0.125001;
}

int main()
{
  struct rv_mesh m = new_box_rv_mesh(3);
  char fname[64];
  for (unsigned it = 0; 1; ++it) {
    struct rv_mesh out = refine_reduced(m, simple);
    printf("%u elements, %u vertices\n", out.nelems, out.nverts);
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
