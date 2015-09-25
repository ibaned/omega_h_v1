#include <math.h>             // for fabs
#include <stdio.h>            // for printf
#include "algebra.h"          // for vector_norm
#include "classify_box.h"     // for mesh_classify_box
#include "coarsen_by_size.h"  // for coarsen_by_size
#include "eval_field.h"       // for mesh_eval_field
#include "mesh.h"             // for mesh_free_nodal_field, free_mesh, mesh_...
#include "quality.h"          // for mesh_min_quality
#include "refine_by_size.h"   // for refine_by_size
#include "vtk.h"              // for write_vtk

static void fine_fun(double const* x, double* s)
{
  double coarse = 0.5;
  double fine = 0.05;
  double radius = vector_norm(x, 3);
  double d = fabs(radius - 0.5);
  s[0] = coarse * d + fine * (1 - d);
}

static void coarse_fun(double const* x, double* s)
{
  (void) x;
  s[0] = 2.0001;
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  char fname[64];
  unsigned it = 0;
  mesh_eval_field(m, 0, "adapt_size", 1, fine_fun);
  while (refine_by_size(&m, 0)) {
    sprintf(fname, "ref_%u.vtu", it++);
    write_vtk(m, fname);
    mesh_free_tag(m, 0, "adapt_size");
    mesh_eval_field(m, 0, "adapt_size", 1, fine_fun);
  }
  mesh_classify_box(m);
  write_vtk(m, "class.vtu");
  double minq = mesh_min_quality(m);
  printf("minq %f\n", minq);
  it = 0;
  mesh_free_tag(m, 0, "adapt_size");
  mesh_eval_field(m, 0, "adapt_size", 1, coarse_fun);
  while (coarsen_by_size(&m, minq, 0.5)) {
    printf("%u elements\n", mesh_count(m, mesh_dim(m)));
    sprintf(fname, "cor_%u.vtu", it++);
    write_vtk(m, fname);
  }
  free_mesh(m);
}
