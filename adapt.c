#include "adapt.h"
#include "refine_by_size.h"
#include "coarsen_by_size.h"
#include "split_slivers.h"
#include "quality.h"
#include "measure_edges.h"
#include "doubles.h"

#include "vtk.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_OPS 30

static void adapt_summary(struct mesh* m)
{
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
  for (unsigned i = 0; i < MAX_OPS; ++i) {
    if (refine_by_size(p_m)) {
      printf("split long edges\n");
      write_vtk_step(*p_m);
      continue;
    }
    if (coarsen_by_size(p_m, mesh_min_quality(*p_m), size_ratio_floor)) {
      printf("collapse short edges\n");
      write_vtk_step(*p_m);
      continue;
    }
    if (mesh_dim(*p_m) < 2) {
      adapt_summary(*p_m);
      return;
    }
    if (split_slivers(p_m, 2, VERT_EDGE_SLIVER, qual_floor, 0)) {
      printf("split vert-edge triangles\n");
      write_vtk_step(*p_m);
      continue;
    }
    if (mesh_dim(*p_m) < 3) {
      adapt_summary(*p_m);
      return;
    }
    if (split_slivers(p_m, 3, EDGE_EDGE_SLIVER, qual_floor, size_ratio_floor)) {
      printf("split edge-edge tets\n");
      write_vtk_step(*p_m);
      continue;
    }
    if (split_slivers(p_m, 3, VERT_FACE_SLIVER, qual_floor, size_ratio_floor)) {
      printf("split vert-face tets\n");
      write_vtk_step(*p_m);
      continue;
    }
  }
  fprintf(stderr, "mesh_adapt could not succeed after %u operations\n", MAX_OPS);
  abort();
}
