#include <math.h>
#include <stdio.h>

#include "algebra.h"
#include "classify_box.h"
#include "coarsen_by_size.h"
#include "derive_model.h"
#include "eval_field.h"
#include "ints.h"
#include "mesh.h"
#include "quality.h"
#include "refine_by_size.h"
#include "subset.h"
#include "vtk.h"

static void fine_fun(double const* x, double* s)
{
  double coarse = 0.5;
  double fine = 0.2;
  double radius = vector_norm(x, 3);
  double d = fabs(radius - 0.5);
  s[0] = coarse * d + fine * (1 - d);
}

static void coarse_fun(double const* x, double* s)
{
  (void) x;
  s[0] = 4;
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  mesh_derive_model(m, PI / 4);
  mesh_set_rep(m, MESH_FULL);
  char fname[64];
  unsigned it = 0;
  mesh_eval_field(m, 0, "adapt_size", 1, fine_fun);
  while (refine_by_size(&m, 0)) {
    sprintf(fname, "ref_%u.vtu", it++);
    write_vtu(m, fname);
    mesh_free_tag(m, 0, "adapt_size");
    mesh_eval_field(m, 0, "adapt_size", 1, fine_fun);
  }
  double minq = 0.1;
  printf("minq %f\n", minq);
  it = 0;
  mesh_free_tag(m, 0, "adapt_size");
  mesh_eval_field(m, 0, "adapt_size", 1, coarse_fun);
  while (coarsen_by_size(&m, minq, 0.5)) {
    printf("%u elements\n", mesh_count(m, mesh_dim(m)));
    sprintf(fname, "cor_%u.vtu", it++);
    write_vtu(m, fname);
  }
  free_mesh(m);
}
