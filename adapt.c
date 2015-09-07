#include "adapt.h"
#include "refine_by_size.h"
#include "coarsen_by_size.h"
#include "refine_slivers.h"
#include "split_slivers.h"
#include "quality.h"
#include "measure_edges.h"
#include "doubles.h"
#include "coarsen_slivers.h"

#include "vtk.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_OPS 30

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

static void satisfy_size(struct mesh** p_m, double size_floor)
{
  double init_qual = mesh_min_quality(*p_m);
  while (refine_by_size(p_m, init_qual))
    incr_op_count(*p_m, "split long edges\n");
  while (coarsen_by_size(p_m, init_qual, size_floor, 0))
    incr_op_count(*p_m, "collapse short edges\n");
}

/*
static struct { unsigned dim; enum sliver_type st; } const sliver_tab[3] = {
  {2, VERT_EDGE_SLIVER},
  {3, EDGE_EDGE_SLIVER},
  {3, VERT_FACE_SLIVER}
};

static char const* const sliver_names[5] = {
  [NOT_SLIVER] = "non",
  [VERT_VERT_SLIVER] = "vert-vert",
  [VERT_EDGE_SLIVER] = "vert-edge",
  [EDGE_EDGE_SLIVER] = "edge-edge",
  [VERT_FACE_SLIVER] = "vert-face"
};
*/

static void satisfy_shape(
    struct mesh** p_m, double size_floor, double qual_floor)
{
  while (1) {
    if (coarsen_slivers(p_m, qual_floor)) {
      incr_op_count(*p_m, "coarsen good sliver verts\n");
      continue;
    }
    if (refine_slivers(p_m, qual_floor, 0, 1)) {
      incr_op_count(*p_m, "split good sliver edges\n");
      continue;
    }
    if (refine_slivers(p_m, qual_floor, 0.1, 0)) {
      incr_op_count(*p_m, "split any sliver edges\n");
      continue;
    }
    if (mesh_min_quality(*p_m) < qual_floor) {
      fprintf(stderr, "ran out of options!\n");
      abort();
    } else {
      return;
    }
  }
}

void mesh_adapt(struct mesh** p_m,
    double size_ratio_floor,
    double qual_floor)
{
  global_op_count = 0;
  adapt_summary(*p_m);
  satisfy_size(p_m, size_ratio_floor);
  satisfy_shape(p_m, size_ratio_floor, qual_floor);
}
