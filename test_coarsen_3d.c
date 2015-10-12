#include <math.h>             // for fabs
#include <stdio.h>            // for printf

#include "algebra.h"          // for vector_norm
#include "classify_box.h"     // for mesh_classify_box
#include "coarsen_by_size.h"  // for coarsen_by_size
#include "eval_field.h"       // for mesh_eval_field
#include "mesh.h"             // for mesh_free_nodal_field, free_mesh, new_b...
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
  s[0] = 2;
}

int main()
{
  struct mesh* m = new_box_mesh(3);
  char fname[64];
  unsigned i = 0;
  mesh_eval_field(m, 0, "adapt_size", 1, fine_fun);
  while (refine_by_size(&m, 0)) {
    sprintf(fname, "ref_%u.vtu", i++);
    write_vtu(m, fname);
    mesh_free_tag(m, 0, "adapt_size");
    mesh_eval_field(m, 0, "adapt_size", 1, fine_fun);
  }
  mesh_classify_box(m);
  write_vtu(m, "init.vtu");
  double init_minq = mesh_min_quality(m);
  printf("init minq %f\n", init_minq);
  double cor_qual_floor = 0.3;
  printf("coarsen quality floor %f\n", cor_qual_floor);
  mesh_free_tag(m, 0, "adapt_size");
  mesh_eval_field(m, 0, "adapt_size", 1, coarse_fun);
  i = 0;
  while (coarsen_by_size(&m, cor_qual_floor, 0.5)) {
    sprintf(fname, "cor_%u.vtu", i++);
    write_vtu(m, fname);
  }
  double cor_minq = mesh_min_quality(m);
  printf("cor minq %f\n", cor_minq);
  free_mesh(m);
}
