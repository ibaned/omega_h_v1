#include "split_slivers.h"
#include "sliver_keys.h"
#include "ints.h"
#include "collect_keys.h"
#include "refine_common.h"
#include <assert.h>
#include <stdlib.h>

/* much debug */
#include "subset.h"
#include "vtk.h"
#include <stdio.h>

unsigned split_slivers(
    struct mesh** p_m,
    unsigned sliver_dim,
    enum sliver_type st,
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
  sliver_keys(sliver_dim, nslivs, verts_of_slivs, coords,
      st, qual_floor, edge_ratio_floor,
      &bad_slivs, &key_of_slivs);
  unsigned something_to_do = ints_max(bad_slivs, nslivs);
  if (!something_to_do) {
    free(bad_slivs);
    free(key_of_slivs);
    return 0;
  }
  unsigned key_dim = 1;
  if (sliver_dim == 3 && st == VERT_FACE_SLIVER)
    key_dim = 2;
  { /* le debug */
    unsigned* sliv_offsets = ints_exscan(bad_slivs, nslivs);
    struct mesh* sm = subset_mesh(m, sliver_dim, sliv_offsets);
    char buf[64];
    static unsigned foo[4][4] = {{}};
    sprintf(buf, "sliv_%u_%u_%u.vtu", sliver_dim, key_dim, foo[sliver_dim][key_dim]++);
    write_vtk(sm, buf);
    free_mesh(sm);
    free(sliv_offsets);
  }
  unsigned nkeys = mesh_count(m, key_dim);
  unsigned const* slivs_of_keys_offsets =
    mesh_ask_up(m, key_dim, sliver_dim)->offsets;
  unsigned const* slivs_of_keys =
    mesh_ask_up(m, key_dim, sliver_dim)->adj;
  unsigned const* slivs_of_keys_directions =
    mesh_ask_up(m, key_dim, sliver_dim)->directions;
  unsigned* candidates = collect_keys(sliver_dim, key_dim, nkeys,
      slivs_of_keys_offsets, slivs_of_keys, slivs_of_keys_directions,
      bad_slivs, key_of_slivs);
  { /* le debug */
    unsigned* key_offsets = ints_exscan(candidates, nkeys);
    struct mesh* sm = subset_mesh(m, key_dim, key_offsets);
    char buf[64];
    static unsigned bar[4][4] = {{}};
    sprintf(buf, "key_%u_%u_%u.vtu", sliver_dim, key_dim, bar[sliver_dim][key_dim]++);
    write_vtk(sm, buf);
    free_mesh(sm);
    free(key_offsets);
  }
  free(bad_slivs);
  free(key_of_slivs);
  refine_common(p_m, key_dim, candidates, 0);
  free(candidates);
  return 1;
}
