#include "reflect_down.h"

#include <assert.h>
#include <stdio.h>

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
   we can use this macro (template) to make it clearer to
   the compiler it should generate unrolled code for
   these two cases.

TODO: consolidate this with find_by_verts.h
*/

#define FIND_LOW_FAST(N) \
LOOP_INOUT static inline unsigned find_low_fast_##N( \
    unsigned const* verts_wanted, \
    unsigned const* verts_of_lows, \
    unsigned const* lows_of_verts_offsets, \
    unsigned const* lows_of_verts) \
{ \
  unsigned v0 = verts_wanted[0]; \
  unsigned f = lows_of_verts_offsets[v0]; \
  unsigned e = lows_of_verts_offsets[v0 + 1]; \
  for (unsigned i = f; i < e; ++i) { \
    unsigned low = lows_of_verts[i]; \
    unsigned const* verts_of_low = verts_of_lows + low * N; \
    unsigned j; \
    for (j = 0; j < N; ++j) { \
      unsigned wanted_vert = verts_wanted[j]; \
      unsigned k; \
      for (k = 0; k < N; ++k) \
        if (verts_of_low[k] == wanted_vert) \
          break; \
      if (k == N) \
        break; \
    } \
    if (j == N) \
      return low; \
  } \
  fprintf(stderr, "could not find new edge (%u %u)\n", \
      verts_wanted[0], verts_wanted[1]); \
  assert(0); \
}

FIND_LOW_FAST(2)
FIND_LOW_FAST(3)

LOOP_INOUT static inline unsigned find_low_fast(
    unsigned nverts_wanted,
    unsigned const* verts_wanted,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts)
{
  if (nverts_wanted == 2)
    return find_low_fast_2(verts_wanted, verts_of_lows, lows_of_verts_offsets,
        lows_of_verts);
  else
    return find_low_fast_3(verts_wanted, verts_of_lows, lows_of_verts_offsets,
        lows_of_verts);
}

LOOP_KERNEL(reflect_down_entity_fast,
    unsigned const* verts_of_highs,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned verts_per_high,
    unsigned lows_per_high,
    unsigned verts_per_low,
    unsigned* lows_of_highs,
    unsigned const* const* high_verts_of_lows)

  unsigned const* verts_of_high = verts_of_highs + i * verts_per_high;
  unsigned* lows_of_high = lows_of_highs + i * lows_per_high;
  for (unsigned j = 0; j < lows_per_high; ++j) {
    unsigned const* high_verts_of_low = high_verts_of_lows[j];
    unsigned verts_wanted[3];
    for (unsigned k = 0; k < verts_per_low; ++k)
      verts_wanted[k] = verts_of_high[high_verts_of_low[k]];
    lows_of_high[j] = find_low_fast(verts_per_low, verts_wanted,
        verts_of_lows, lows_of_verts_offsets, lows_of_verts);
  }
}

static unsigned* reflect_down_fast(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* verts_of_highs,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts)
{
  unsigned verts_per_high = the_down_degrees[high_dim][0];
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned verts_per_low = the_down_degrees[low_dim][0];
  assert(verts_per_low == 2 || verts_per_low == 3);
  unsigned* lows_of_highs = LOOP_MALLOC(unsigned, nhighs * lows_per_high);
  unsigned const* const* high_verts_of_lows =
    the_canonical_orders[high_dim][low_dim][0];
  LOOP_EXEC(reflect_down_entity_fast, nhighs,
      verts_of_highs,
      verts_of_lows,
      lows_of_verts_offsets,
      lows_of_verts,
      verts_per_high,
      lows_per_high,
      verts_per_low,
      lows_of_highs,
      high_verts_of_lows);
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
  return reflect_down_fast(high_dim, low_dim, nhighs,
      verts_of_highs, verts_of_lows, lows_of_verts->offsets, lows_of_verts->adj);
}
