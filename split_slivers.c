#include "split_slivers.h"
#include "bad_elem_keys.h"
#include "ints.h"
#include "collect_keys.h"
#include "refine_common.h"
#include <assert.h>
#include <stdlib.h>

unsigned split_slivers(
    unsigned sliver_dim,
    struct mesh** p_m,
    double qual_floor,
    double edge_ratio_floor)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  assert(elem_dim >= 2);
  assert(elem_dim <= 3);
  unsigned const* verts_of_slivs = mesh_ask_down(m, sliver_dim, 0);
  unsigned nslivs = mesh_count(m, sliver_dim);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  unsigned* bad_slivs;
  unsigned* key_of_slivs;
  bad_elem_keys(2, nslivs, verts_of_slivs, coords,
      SLIVER_ELEM, qual_floor, edge_ratio_floor,
      &bad_slivs, &key_of_slivs);
  unsigned something_to_do = ints_max(bad_slivs, nslivs);
  if (!something_to_do) {
    free(bad_slivs);
    free(key_of_slivs);
    return 0;
  }
  unsigned nedges = mesh_count(m, 1);
  unsigned const* slivs_of_edges_offsets = mesh_ask_up(m, 1, sliver_dim)->offsets;
  unsigned const* slivs_of_edges = mesh_ask_up(m, 1, sliver_dim)->adj;
  unsigned const* slivs_of_edges_directions = mesh_ask_up(m, 1, sliver_dim)->directions;
  unsigned* candidates = collect_keys(sliver_dim, 1, nedges,
      slivs_of_edges_offsets, slivs_of_edges, slivs_of_edges_directions,
      bad_slivs, key_of_slivs);
  free(bad_slivs);
  free(key_of_slivs);
  refine_common(p_m, 1, candidates);
  free(candidates);
  return 1;
}
