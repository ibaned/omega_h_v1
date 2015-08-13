#include "tables.h"
#include "refine_by_size.h"
#include "classif_box.h"
#include "coarsen_by_size.h"
#include "vtk.h"
#include "verify.h"
#include "algebra.h"
#include "element_qualities.h"
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
  unsigned elem_dim = 2;
  unsigned nelems;
  unsigned nverts;
  unsigned* verts_of_elems;
  double* coords;
  get_box_copy(elem_dim, &nelems, &nverts, &verts_of_elems, &coords);
  while (refine_by_size(elem_dim, &nelems, &nverts,
             &verts_of_elems, &coords, 0, fine_fun));
  unsigned* class_dim = classif_box(elem_dim, nverts, coords);
  write_vtk("init.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
//double minq = min_element_quality(elem_dim, nelems, verts_of_elems, coords);
  double minq = 0.001;
  printf("minq %f\n", minq);
  coarsen_by_size(elem_dim, &nelems, &nverts,
      &verts_of_elems, &coords, &class_dim, coarse_fun, minq);
  coarsen_by_size(elem_dim, &nelems, &nverts,
      &verts_of_elems, &coords, &class_dim, coarse_fun, minq);
  write_vtk("pre_cor.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
  minq = min_element_quality(elem_dim, nelems, verts_of_elems, coords);
  printf("cor minq %f\n", minq);
  double qual_floor = 0.5;
  double edge_ratio_floor = 0; /* never rule "short edge" */
  split_sliver_tris(elem_dim, &nelems, &nverts, &verts_of_elems, &coords,
      &class_dim, qual_floor, edge_ratio_floor);
  write_vtk("sliver.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
  coarsen_by_size(elem_dim, &nelems, &nverts,
      &verts_of_elems, &coords, &class_dim, coarse_fun, minq);
  write_vtk("post_cor.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
  free(verts_of_elems);
  free(coords);
  free(class_dim);
}
