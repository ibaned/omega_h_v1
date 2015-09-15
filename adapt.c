#include "adapt.h"
#include <stdio.h>            // for printf, fprintf, stderr
#include <stdlib.h>           // for abort
#include "loop.h"             // for loop_free
#include "coarsen_by_size.h"  // for coarsen_by_size
#include "coarsen_slivers.h"  // for coarsen_slivers
#include "doubles.h"          // for doubles_max, doubles_min
#include "field.h"            // for const_field
#include "measure_edges.h"    // for measure_edges
#include "mesh.h"             // for mesh_count, mesh_dim, mesh_find_nodal_f...
#include "quality.h"          // for mesh_min_quality
#include "refine_by_size.h"   // for refine_by_size
#include "size.h"             // for mesh_domain_size
#include "swap_slivers.h"     // for swap_slivers
#include "vtk.h"              // for write_vtk_step

static unsigned global_op_count = 0;
static unsigned global_max_ops = 0;

static void adapt_summary(struct mesh* m)
{
  printf("%u elements, ", mesh_count(m, mesh_dim(m)));
  printf("min quality %.2e, ", mesh_min_quality(m));
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  double const* size = mesh_find_nodal_field(m, "adapt_size")->data;
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, size);
  double min = doubles_min(edge_sizes, nedges);
  double max = doubles_max(edge_sizes, nedges);
  loop_free(edge_sizes);
  printf("metric range %.2e - %.2e ", max, min);
  printf("domain size %.6e\n", mesh_domain_size(m));
}

static void incr_op_count(struct mesh* m, char const* what)
{
  if (global_op_count > global_max_ops) {
    fprintf(stderr, "mesh_adapt could not succeed after %u operations\n",
        global_max_ops);
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
  while (coarsen_by_size(p_m, qual_floor, size_floor))
    incr_op_count(*p_m, "collapse short edges\n");
}

static void satisfy_shape(
    struct mesh** p_m,
    double qual_floor,
    unsigned nsliver_layers)
{
  while (1) {
    double prev_qual = mesh_min_quality(*p_m);
    if (prev_qual >= qual_floor)
      return;
    if (mesh_dim(*p_m) == 3 &&
        swap_slivers(p_m, qual_floor, nsliver_layers)) {
      incr_op_count(*p_m, "swap good edges\n");
      continue;
    }
    if (coarsen_slivers(p_m, qual_floor, nsliver_layers)) {
      incr_op_count(*p_m, "coarsen good verts\n");
      continue;
    }
    fprintf(stderr, "ran out of options!\n");
    abort();
  }
}

void mesh_adapt(struct mesh** p_m,
    double size_ratio_floor,
    double good_qual,
    unsigned nsliver_layers,
    unsigned max_ops)
{
  global_op_count = 0;
  global_max_ops = max_ops;
  adapt_summary(*p_m);
  satisfy_size(p_m, size_ratio_floor, good_qual);
  satisfy_shape(p_m, good_qual, nsliver_layers);
}
