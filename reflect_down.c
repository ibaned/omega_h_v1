#include "reflect_down.h"

#include <assert.h>

#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

/* some optimization here due to the relatively
   high cost of these operations.
   first, nverts_wanted is either 2 or 3.
   that is because reflect_down is only
   used to derive intermediate downward adjacencies,
   and we always have the n->0 adjacency, so it must
   be one of:
     3->2
     3->1
     2->1
   and for d=2, d=1, the number of vertices of a simplex
   is n=3, n=2, respectively.

TODO: consolidate this with find_by_verts.h
*/

LOOP_INOUT static inline unsigned find_face_up(
    unsigned const* verts_wanted,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned const* lows_of_verts_directions)
{
  unsigned v0 = verts_wanted[0];
  unsigned f = lows_of_verts_offsets[v0];
  unsigned e = lows_of_verts_offsets[v0 + 1];
  for (unsigned i = f; i < e; ++i) {
    unsigned low = lows_of_verts[i];
    unsigned dir = lows_of_verts_directions[i];
    unsigned const* verts_of_low = verts_of_lows + low * 3;
    if (((verts_of_low[(dir + 1) % 3] == verts_wanted[1]) &&
         (verts_of_low[(dir + 2) % 3] == verts_wanted[2])) ||
        ((verts_of_low[(dir + 1) % 3] == verts_wanted[2]) &&
         (verts_of_low[(dir + 2) % 3] == verts_wanted[1])))
      return low;
  }
  LOOP_NORETURN(INVALID);
}

LOOP_INOUT static inline unsigned find_edge_up(
    unsigned const* verts_wanted,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned const* lows_of_verts_directions)
{
  unsigned v0 = verts_wanted[0];
  unsigned f = lows_of_verts_offsets[v0];
  unsigned e = lows_of_verts_offsets[v0 + 1];
  for (unsigned i = f; i < e; ++i) {
    unsigned low = lows_of_verts[i];
    unsigned dir = lows_of_verts_directions[i];
    unsigned const* verts_of_low = verts_of_lows + low * 2;
    if (verts_of_low[1 - dir] == verts_wanted[1])
      return low;
  }
  LOOP_NORETURN(INVALID);
}

LOOP_INOUT static inline unsigned find_low(
    unsigned nverts_wanted,
    unsigned const* verts_wanted,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned const* lows_of_verts_directions)
{
  if (nverts_wanted == 2)
    return find_edge_up(verts_wanted, verts_of_lows,
        lows_of_verts_offsets, lows_of_verts, lows_of_verts_directions);
  else
    return find_face_up(verts_wanted, verts_of_lows,
        lows_of_verts_offsets, lows_of_verts, lows_of_verts_directions);
}

LOOP_KERNEL(reflect_down_entity,
    unsigned const* verts_of_highs,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned const* lows_of_verts_directions,
    unsigned verts_per_high,
    unsigned lows_per_high,
    unsigned verts_per_low,
    unsigned* lows_of_highs,
    unsigned** high_verts_of_lows)

  unsigned const* verts_of_high = verts_of_highs + i * verts_per_high;
  unsigned* lows_of_high = lows_of_highs + i * lows_per_high;
  for (unsigned j = 0; j < lows_per_high; ++j) {
    unsigned const* high_verts_of_low = high_verts_of_lows[j];
    unsigned verts_wanted[3];
    for (unsigned k = 0; k < verts_per_low; ++k)
      verts_wanted[k] = verts_of_high[high_verts_of_low[k]];
    lows_of_high[j] = find_low(verts_per_low, verts_wanted, verts_of_lows,
        lows_of_verts_offsets, lows_of_verts, lows_of_verts_directions);
  }
}

static unsigned* reflect_down(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* verts_of_highs,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned const* lows_of_verts_directions)
{
  unsigned verts_per_high = the_down_degrees[high_dim][0];
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned verts_per_low = the_down_degrees[low_dim][0];
  assert(verts_per_low == 2 || verts_per_low == 3);
  unsigned* lows_of_highs = LOOP_MALLOC(unsigned, nhighs * lows_per_high);
  unsigned** high_verts_of_lows = orders_to_device(high_dim, low_dim, 0);
  LOOP_EXEC(reflect_down_entity, nhighs,
      verts_of_highs,
      verts_of_lows,
      lows_of_verts_offsets,
      lows_of_verts,
      lows_of_verts_directions,
      verts_per_high,
      lows_per_high,
      verts_per_low,
      lows_of_highs,
      high_verts_of_lows);
  free_orders(high_verts_of_lows, high_dim, low_dim);
  return lows_of_highs;
}

unsigned* mesh_reflect_down(
    struct mesh* m,
    unsigned high_dim,
    unsigned low_dim)
{
  unsigned nhighs = mesh_count(m, high_dim);
  unsigned const* verts_of_highs = mesh_ask_down(m, high_dim, 0);
  struct const_up* lows_of_verts = mesh_ask_up(m, 0, low_dim);
  unsigned const* verts_of_lows = mesh_ask_down(m, low_dim, 0);
  return reflect_down(high_dim, low_dim, nhighs, verts_of_highs, verts_of_lows,
      lows_of_verts->offsets, lows_of_verts->adj, lows_of_verts->directions);
}
