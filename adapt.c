#include "adapt.h"
#include "refine_by_size.h"
#include "coarsen_by_size.h"
#include "split_slivers.h"

#include "vtk.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_OPS 30

void mesh_adapt(struct mesh** p_m,
    double size_ratio_floor,
    double qual_floor)
{
  for (unsigned i = 0; i < MAX_OPS; ++i) {
    if (refine_by_size(p_m)) {
      printf("refine long edges\n");
      write_vtk_step(*p_m);
      continue;
    }
    if (coarsen_by_size(p_m, mesh_min_quality(*p_m), size_ratio_floor)) {
      printf("coarsen edges < %f h\n", size_ratio_floor);
      write_vtk_step(*p_m);
      continue;
    }
    if (mesh_dim(*p_m) < 2)
      return;
    if (split_slivers(p_m, 2, VERT_EDGE_SLIVER, qual_floor, size_ratio_floor)) {
      printf("split vert-edge triangles\n");
      write_vtk_step(*p_m);
      continue;
    }
    if (mesh_dim(*p_m) < 3)
      return;
    if (split_slivers(p_m, 3, VERT_EDGE_SLIVER, qual_floor, size_ratio_floor)) {
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
