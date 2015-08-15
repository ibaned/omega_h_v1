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
  unsigned elem_dim = 3;
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
  double init_minq = min_element_quality(elem_dim, nelems, verts_of_elems, coords);
  printf("init minq %f\n", init_minq);
  double cor_qual_floor = 0.4;
  printf("coarsen quality floor %f\n", cor_qual_floor);
  while (coarsen_by_size(elem_dim, &nelems, &nverts,
      &verts_of_elems, &coords, &class_dim, coarse_fun, cor_qual_floor));
  write_vtk("cor.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
  double cor_minq = min_element_quality(elem_dim, nelems, verts_of_elems, coords);
  printf("cor minq %f\n", cor_minq);
  double sliver_qual_floor = 0.5;
  printf("sliver quality floor %f\n", sliver_qual_floor);
  double edge_ratio_floor = 0; /* never rule "short edge" */
  split_sliver_tris(elem_dim, &nelems, &nverts, &verts_of_elems, &coords,
      &class_dim, sliver_qual_floor, edge_ratio_floor);
  write_vtk("sliver.vtu", elem_dim, nelems, nverts, verts_of_elems, coords,
      class_dim);
  double sliver_minq = min_element_quality(elem_dim, nelems, verts_of_elems, coords);
  printf("sliver minq %f\n", sliver_minq);
  free(verts_of_elems);
  free(coords);
  free(class_dim);
}
