#include "warp_to_limit.h"
#include "doubles.h"
#include "quality.h"
#include "mesh.h"
#include <stdlib.h>
#include <stdio.h>

unsigned warp_to_limit(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords,
    double const* warps,
    double** p_coords,
    double** p_warps)
{
  unsigned hit_limit = 0;
  double* coords_out = malloc(sizeof(double) * 3 * nverts);
  double factor = 1;
  double* warps_out = 0;
  doubles_axpy(factor, warps, coords, coords_out, 3 * nverts);
  while (min_element_quality(
        elem_dim, nelems, verts_of_elems, coords_out) <= 0) {
    hit_limit = 1;
    factor /= 2;
    doubles_axpy(factor, warps, coords, coords_out, 3 * nverts);
  }
  *p_coords = coords_out;
  if (p_warps) {
    warps_out = malloc(sizeof(double) * 3 * nverts);
    doubles_axpy(-factor, warps, warps, warps_out, 3 * nverts);
    *p_warps = warps_out;
  }
  return !hit_limit;
}

unsigned mesh_warp_to_limit(struct mesh* m)
{
  double* coords;
  double* warps;
  unsigned ok =
      warp_to_limit(mesh_dim(m), mesh_count(m, mesh_dim(m)), mesh_count(m, 0),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_nodal_field(m, "coordinates")->data,
      mesh_find_nodal_field(m, "warp")->data,
      &coords, &warps);
  mesh_free_nodal_field(m, "coordinates");
  mesh_free_nodal_field(m, "warp");
  mesh_add_nodal_field(m, "coordinates", 3, coords);
  mesh_add_nodal_field(m, "warp", 3, warps);
  return ok;
}
