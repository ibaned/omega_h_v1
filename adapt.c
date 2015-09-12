#include "adapt.h"
#include "mesh.h"
#include "refine_by_size.h"
#include "coarsen_by_size.h"
#include "refine_slivers.h"
#include "split_slivers.h"
#include "quality.h"
#include "measure_edges.h"
#include "doubles.h"
#include "coarsen_slivers.h"
#include "swap_slivers.h"

#include "vtk.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_OPS 100

static unsigned global_op_count = 0;

static void adapt_summary(struct mesh* m)
{
  printf("%u elements, ", mesh_count(m, mesh_dim(m)));
  printf("min quality %f, ", mesh_min_quality(m));
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  double const* size = mesh_find_nodal_field(m, "adapt_size")->data;
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, size);
  double min = doubles_min(edge_sizes, nedges);
  double max = doubles_max(edge_sizes, nedges);
  free(edge_sizes);
  printf("metric range %f - %f\n", max, min);
}

static void incr_op_count(struct mesh* m, char const* what)
{
  if (global_op_count > MAX_OPS) {
    fprintf(stderr, "mesh_adapt could not succeed after %u operations\n", MAX_OPS);
    abort();
  }
  ++global_op_count;
  printf("%s", what);
  adapt_summary(m);
  write_vtk_step(m);
}

static void satisfy_size(struct mesh** p_m, double size_floor, double good_qual)
{
  double qual_floor = mesh_min_quality(*p_m);
  if (good_qual < qual_floor)
    qual_floor = good_qual;
  while (refine_by_size(p_m, qual_floor))
    incr_op_count(*p_m, "split long edges\n");
  while (coarsen_by_size(p_m, qual_floor, size_floor, 0))
    incr_op_count(*p_m, "collapse short edges\n");
}

static void satisfy_shape(
    struct mesh** p_m, double qual_floor)
{
  while (1) {
    double prev_qual = mesh_min_quality(*p_m);
    if (prev_qual >= qual_floor)
      return;
    if (mesh_dim(*p_m) == 3 && swap_slivers(p_m, 1.0, 0.0, 0)) {
      incr_op_count(*p_m, "swap good edges\n");
      abort();
      continue;
    }
    if (coarsen_slivers(p_m, 1.0, 0)) {
      incr_op_count(*p_m, "coarsen good verts\n");
      continue;
    }
    if (refine_slivers(p_m, 1, 1.0, 0.0, 1, 0)) {
      incr_op_count(*p_m, "split good edges\n");
      continue;
    }
    fprintf(stderr, "ran out of options!\n");
    abort();
  }
}

void mesh_adapt(struct mesh** p_m,
    double size_ratio_floor,
    double good_qual)
{
  global_op_count = 0;
  adapt_summary(*p_m);
  satisfy_size(p_m, size_ratio_floor, good_qual);
  satisfy_shape(p_m, good_qual);
}
