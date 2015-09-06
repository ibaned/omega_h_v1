#include "adapt.h"
#include "refine_by_size.h"
#include "coarsen_by_size.h"
#include "refine_slivers.h"
#include "quality.h"
#include "measure_edges.h"
#include "doubles.h"

#include "vtk.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_OPS 50

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

void mesh_adapt(struct mesh** p_m,
    double size_ratio_floor,
    double qual_floor)
{
  double prev_qual = mesh_min_quality(*p_m);
  for (unsigned i = 0; i < MAX_OPS; ++i) {
    if (refine_by_size(p_m)) {
      printf("split long edges\n");
      write_vtk_step(*p_m);
      prev_qual = mesh_min_quality(*p_m);
      continue;
    }
    double coarsen_qual = prev_qual + 1e-10;
    if (qual_floor < prev_qual)
      coarsen_qual = qual_floor;
    printf("coarsen qual floor %f\n", coarsen_qual);
    if (coarsen_by_size(p_m, coarsen_qual, size_ratio_floor)) {
      printf("collapse short edges\n");
      write_vtk_step(*p_m);
      prev_qual = mesh_min_quality(*p_m);
      printf("qual after coarsen %f\n", prev_qual);
      continue;
    }
    if (refine_slivers(p_m, qual_floor)) {
      printf("split sliver edges\n");
      write_vtk_step(*p_m);
      printf("qual after sliver split %f\n", mesh_min_quality(*p_m));
      continue;
    }
    adapt_summary(*p_m);
    return;
  }
  fprintf(stderr, "mesh_adapt could not succeed after %u operations\n", MAX_OPS);
  abort();
}
