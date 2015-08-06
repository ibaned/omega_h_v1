#include "refine_reduced.h"
#include "tables.h"
#include "vtk.h"
#include <stdlib.h>
#include <stdio.h>

static double sf(double const x[])
{
  (void) x;
  return 1.0;
}

int main()
{
  struct rv_mesh in = new_box_rv_mesh(3);
  struct rv_mesh out = refine_reduced(in, sf);
  free_rv_mesh(in);
  printf("%u elements, %u vertices\n", out.nelems, out.nverts);
  write_vtk("refined.vtu", out.elem_dim, out.nelems, out.nverts,
      out.verts_of_elems, out.coords);
  free_rv_mesh(out);
}
