#include "tables.h"
#include "refine_reduced.h"
#include "classif_box.h"
#include "coarsen_reduced.h"
#include "vtk.h"
#include "verify.h"
#include "algebra.h"
#include <stdio.h>
#include <stdlib.h>

static double fine_fun(double const x[])
{
  double coarse = 0.5;
  double fine = 0.05;
  double radius = vector_norm(x, 3);
  double d = fabs(radius - 0.5);
  return coarse * d + fine * (1 - d);
}

static double coarse_fun(double const x[])
{
  (void) x;
  return 2.0001;
}

int main()
{
  unsigned elem_dim = 3;
  unsigned nelems;
  unsigned nverts;
  unsigned* verts_of_elems;
  double* coords;
  get_box_copy(elem_dim, &nelems, &nverts, &verts_of_elems, &coords);
  char fname[64];
  unsigned it = 0;
  while (refine_reduced(elem_dim, &nelems, &nverts,
             &verts_of_elems, &coords, fine_fun)) {
    sprintf(fname, "ref_%u.vtu", it++);
    write_vtk(fname, elem_dim, nelems, nverts, verts_of_elems, coords, 0);
  }
  unsigned* class_dim = classif_box(elem_dim, nverts, coords);
  write_vtk("class.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
  it = 0;
  while (coarsen_reduced(elem_dim, &nelems, &nverts,
             &verts_of_elems, &coords, &class_dim, coarse_fun)) {
    verify(elem_dim, nelems, nverts, verts_of_elems, coords);
    sprintf(fname, "cor_%u.vtu", it++);
    write_vtk(fname, elem_dim, nelems, nverts, verts_of_elems, coords,
        class_dim);
  }
  free(verts_of_elems);
  free(coords);
  free(class_dim);
}
