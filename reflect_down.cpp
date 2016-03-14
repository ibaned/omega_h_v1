#include "reflect_down.hpp"

#include <cassert>

#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tables.hpp"
#include "find_by_verts.hpp"

namespace omega_h {

LOOP_KERNEL(reflect_down_entity,
    unsigned high_dim,
    unsigned low_dim,
    unsigned const* verts_of_highs,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned const* lows_of_verts_directions,
    unsigned verts_per_high,
    unsigned lows_per_high,
    unsigned verts_per_low,
    unsigned* lows_of_highs)

  unsigned const* const* high_verts_of_lows =
      the_canonical_orders[high_dim][low_dim][0];
  unsigned const* verts_of_high = verts_of_highs + i * verts_per_high;
  unsigned* lows_of_high = lows_of_highs + i * lows_per_high;
  for (unsigned j = 0; j < lows_per_high; ++j) {
    unsigned const* high_verts_of_low = high_verts_of_lows[j];
    unsigned verts_wanted[3];
    for (unsigned k = 0; k < verts_per_low; ++k)
      verts_wanted[k] = verts_of_high[high_verts_of_low[k]];
    lows_of_high[j] = find_by_verts(low_dim, verts_wanted, verts_of_lows,
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
  LOOP_EXEC(reflect_down_entity, nhighs,
      high_dim,
      low_dim,
      verts_of_highs,
      verts_of_lows,
      lows_of_verts_offsets,
      lows_of_verts,
      lows_of_verts_directions,
      verts_per_high,
      lows_per_high,
      verts_per_low,
      lows_of_highs);
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

}
