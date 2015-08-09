#include "tables.h"
#include "refine_reduced.h"
#include "classif_box.h"
#include "coarsen_reduced.h"
#include "vtk.h"
#include <stdio.h>
#include <stdlib.h>

static double fine(double const x[])
{
  (void) x;
  return 0.5001;
}

static double coarse(double const x[])
{
  (void) x;
  return 2.0001;
}

int main()
{
  unsigned elem_dim = 2;
  unsigned nelems;
  unsigned nverts;
  unsigned* verts_of_elems;
  double* coords;
  get_box_copy(elem_dim, &nelems, &nverts, &verts_of_elems, &coords);
  char fname[64];
  unsigned it = 0;
  while (refine_reduced(elem_dim, &nelems, &nverts,
             &verts_of_elems, &coords, fine)) {
    sprintf(fname, "ref_%u.vtu", it++);
    write_vtk(fname, elem_dim, nelems, nverts, verts_of_elems, coords, 0);
  }
  unsigned* class_dim = classif_box(elem_dim, nverts, coords);
  write_vtk("class.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
  coarsen_reduced(elem_dim, &nelems, &nverts, &verts_of_elems,
      &coords, &class_dim, coarse);
  write_vtk("coarse.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
  free(verts_of_elems);
  free(coords);
  free(class_dim);
}
