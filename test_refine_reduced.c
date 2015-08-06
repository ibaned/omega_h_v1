#include "refine_reduced.h"
#include "tables.h"
#include "vtk.h"
#include <stdlib.h>
#include <stdio.h>

static double sf(double const x[])
{
  (void) x;
  return 0.501;
}

int main()
{
  struct rv_mesh m = new_box_rv_mesh(3);
  while (1) {
    struct rv_mesh out = refine_reduced(m, sf);
    printf("%u elements, %u vertices\n", out.nelems, out.nverts);
    if (out.nelems == m.nelems) {
      free_rv_mesh(out);
      break;
    }
    free_rv_mesh(m);
    m = out;
  }
  write_vtk("refined.vtu", m.elem_dim, m.nelems, m.nverts,
      m.verts_of_elems, m.coords);
  free_rv_mesh(m);
}
