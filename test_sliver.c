#include "tables.h"
#include "refine_by_size.h"
#include "classify_box.h"
#include "coarsen_by_size.h"
#include "vtk.h"
#include "verify.h"
#include "algebra.h"
#include "quality.h"
#include "split_sliver_tris.h"
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
  return 2;
}

int main()
{
  struct mesh* m = new_box_mesh(3);
  char fname[64];
  unsigned i = 0;
  while (refine_by_size(&m, fine_fun)) {
    sprintf(fname, "ref_%u.vtu", i++);
    write_vtk(m, fname);
  }
  mesh_classify_box(m);
  write_vtk(m, "init.vtu");
  double init_minq = mesh_min_quality(m);
  printf("init minq %f\n", init_minq);
  double cor_qual_floor = 0.3;
  printf("coarsen quality floor %f\n", cor_qual_floor);
  i = 0;
  while (coarsen_by_size(&m, coarse_fun, cor_qual_floor)) {
    sprintf(fname, "cor_%u.vtu", i++);
    write_vtk(m, fname);
  }
  double cor_minq = mesh_min_quality(m);
  printf("cor minq %f\n", cor_minq);
  free_mesh(m);
}
