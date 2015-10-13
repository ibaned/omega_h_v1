#include "warp_to_limit.h"

#include <stdio.h>

#include "doubles.h"
#include "loop.h"
#include "mesh.h"
#include "quality.h"
#include "tag.h"

static unsigned exceeds_limit(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords,
    double qual_floor)
{
  double qual = min_element_quality(elem_dim, nelems,
      verts_of_elems, coords);
  if (qual < qual_floor) {
    printf("quality %f exceeds limit\n", qual);
    return 1;
  }
  return 0;
}

unsigned warp_to_limit(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords,
    double const* warps,
    double qual_floor,
    double** p_coords,
    double** p_warps)
{
  unsigned hit_limit = 0;
  double* coords_out = LOOP_MALLOC(double, 3 * nverts);
  double factor = 1;
  double* warps_out = 0;
  doubles_axpy(factor, warps, coords, coords_out, 3 * nverts);
  while (exceeds_limit(elem_dim, nelems,
        verts_of_elems, coords_out, qual_floor)) {
    hit_limit = 1;
    factor /= 2;
    doubles_axpy(factor, warps, coords, coords_out, 3 * nverts);
  }
  *p_coords = coords_out;
  if (p_warps) {
    warps_out = LOOP_MALLOC(double, 3 * nverts);
    doubles_axpy(-factor, warps, warps, warps_out, 3 * nverts);
    *p_warps = warps_out;
  }
  return !hit_limit;
}

unsigned mesh_warp_to_limit(struct mesh* m, double qual_floor)
{
  double* coords;
  double* warps;
  unsigned ok = warp_to_limit(
      mesh_dim(m),
      mesh_count(m, mesh_dim(m)),
      mesh_count(m, 0),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_tag(m, 0, "coordinates")->data,
      mesh_find_tag(m, 0, "warp")->data,
      qual_floor,
      &coords, &warps);
  mesh_free_tag(m, 0, "coordinates");
  mesh_free_tag(m, 0, "warp");
  mesh_add_tag(m, 0, TAG_F64, "coordinates", 3, coords);
  mesh_add_tag(m, 0, TAG_F64, "warp", 3, warps);
  return ok;
}
