#include "tables.h"
#include "refine_by_size.h"
#include "classify_box.h"
#include "coarsen_by_size.h"
#include "vtk.h"
#include "algebra.h"
#include "quality.h"
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
  struct mesh* m = new_box_mesh(2);
  char fname[64];
  unsigned it = 0;
  while (refine_by_size(&m, fine_fun)) {
    sprintf(fname, "ref_%u.vtu", it++);
    write_vtk(m, fname);
  }
  mesh_classify_box(m);
  write_vtk(m, "class.vtu");
  double minq = mesh_min_quality(m);
  printf("minq %f\n", minq);
  it = 0;
  while (coarsen_by_size(&m, coarse_fun, minq, 0.5)) {
    printf("%u elements\n", mesh_count(m, mesh_dim(m)));
    sprintf(fname, "cor_%u.vtu", it++);
    write_vtk(m, fname);
  }
  free_mesh(m);
}
