#include "split_sliver_tris.h"
#include "bad_elem_keys.h"
#include "ints.h"
#include "collect_keys.h"
#include "refine_common.h"
#include <assert.h>
#include <stdlib.h>

unsigned split_sliver_tris(
    struct mesh** p_m,
    double qual_floor,
    double edge_ratio_floor)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  assert(elem_dim >= 2);
  assert(elem_dim <= 3);
  unsigned const* verts_of_tris = mesh_ask_down(m, 2, 0);
  unsigned ntris = mesh_count(m, 2);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  unsigned* bad_tris;
  unsigned* key_of_tris;
  bad_elem_keys(2, ntris, verts_of_tris, coords,
      SLIVER_ELEM, qual_floor, edge_ratio_floor,
      &bad_tris, &key_of_tris);
  unsigned something_to_do = ints_max(bad_tris, ntris);
  if (!something_to_do) {
    free(bad_tris);
    free(key_of_tris);
    return 0;
  }
  unsigned nedges = mesh_count(m, 1);
  unsigned const* tris_of_edges_offsets = mesh_ask_up(m, 1, 2)->offsets;
  unsigned const* tris_of_edges = mesh_ask_up(m, 1, 2)->adj;
  unsigned const* tris_of_edges_directions = mesh_ask_up(m, 1, 2)->directions;
  unsigned* candidates = collect_keys(2, 1, nedges,
      tris_of_edges_offsets, tris_of_edges, tris_of_edges_directions,
      bad_tris, key_of_tris);
  free(bad_tris);
  free(key_of_tris);
  refine_common(p_m, 1, candidates);
  free(candidates);
  return 1;
}
