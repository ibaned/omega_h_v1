#include "adapt.h"
#include "refine_by_size.h"
#include "coarsen_by_size.h"
#include "refine_slivers.h"
#include "split_slivers.h"
#include "quality.h"
#include "measure_edges.h"
#include "doubles.h"

#include "vtk.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_OPS 30

static unsigned global_op_count = 0;

static void incr_op_count(void)
{
  if (global_op_count > MAX_OPS) {
    fprintf(stderr, "mesh_adapt could not succeed after %u operations\n", MAX_OPS);
    abort();
  }
  ++global_op_count;
}

static void adapt_summary(struct mesh* m)
{
  printf("%u elements\n", mesh_count(m, mesh_dim(m)));
  printf("minimum element quality %f\n", mesh_min_quality(m));
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  double const* size = mesh_find_nodal_field(m, "adapt_size")->data;
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, size);
  double min = doubles_min(edge_sizes, nedges);
  double max = doubles_max(edge_sizes, nedges);
  free(edge_sizes);
  printf("edge metric range %f - %f\n", max, min);
}

static void satisfy_size(struct mesh** p_m, double size_floor)
{
  double init_qual = mesh_min_quality(*p_m);
  while (refine_by_size(p_m, init_qual)) {
    incr_op_count();
    printf("split long edges\n");
    write_vtk_step(*p_m);
  }
  while (coarsen_by_size(p_m, init_qual, size_floor)) {
    incr_op_count();
    printf("collapse short edges\n");
    write_vtk_step(*p_m);
  }
}

static void satisfy_shape(struct mesh** p_m, double size_floor, double qual_floor)
{
  while (1) {
    double prev_qual = mesh_min_quality(*p_m);
    if (prev_qual >= qual_floor)
      return;
    if (!split_slivers(p_m, 2, VERT_EDGE_SLIVER, qual_floor, 0)) {
      fprintf(stderr, "couldn't apply sliver splits!\n");
      abort();
    } else {
      incr_op_count();
      printf("split sliver edges\n");
      write_vtk_step(*p_m);
    }
    if (coarsen_by_size(p_m, prev_qual + 1e-10, size_floor)) {
      incr_op_count();
      printf("collapse short (sliver) edges\n");
      write_vtk_step(*p_m);
    }
  }
}

void mesh_adapt(struct mesh** p_m,
    double size_ratio_floor,
    double qual_floor)
{
  global_op_count = 0;
  satisfy_size(p_m, size_ratio_floor);
  satisfy_shape(p_m, size_ratio_floor, qual_floor);
  adapt_summary(*p_m);
}
