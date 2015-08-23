#include "tables.h"
#include "refine_by_size.h"
#include "classify_box.h"
#include "coarsen_by_size.h"
#include "vtk.h"
#include "verify.h"
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
  mesh_add_nodal_label(m, "class_dim",
      classify_box(mesh_dim(m), mesh_count(m, 0),
        mesh_find_nodal_field(m, "coordinates")->data));
  write_vtk(m, "class.vtu");
  double minq = min_element_quality(mesh_dim(m),
      mesh_count(m, mesh_dim(m)),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_nodal_field(m, "coordinates")->data);
  printf("minq %f\n", minq);
  it = 0;
  while (coarsen_by_size(&m, coarse_fun, minq)) {
    verify(mesh_dim(m), mesh_count(m, mesh_dim(m)),
        mesh_count(m, 0),
        mesh_ask_down(m, mesh_dim(m), 0),
        mesh_find_nodal_field(m, "coordinates")->data);
    printf("%u elements\n", mesh_count(m, mesh_dim(m)));
    sprintf(fname, "cor_%u.vtu", it++);
    write_vtk(m, fname);
  }
  free_mesh(m);
}
