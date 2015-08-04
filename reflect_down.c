#include "reflect_down.h"
#include "intersect.h"
#include "tables.h"
#include <stdlib.h>
#include <assert.h>

unsigned* reflect_down(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned nlows,
    unsigned const* verts_of_highs,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts)
{
  unsigned verts_per_high = the_down_degrees[high_dim][0];
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned verts_per_low = the_down_degrees[low_dim][0];
  unsigned* lows_of_highs = malloc(sizeof(unsigned) * nlows * verts_per_low);
  unsigned const* const* high_verts_of_lows =
    the_canonical_orders[high_dim][low_dim][0];
  for (unsigned i = 0; i < nhighs; ++i) {
    unsigned const* up_vert = verts_of_highs + i * verts_per_high;
    unsigned* lows_of_high = lows_of_highs + i * lows_per_high;
    for (unsigned j = 0; j < lows_per_high; ++j) {
      unsigned const* high_verts_of_low = high_verts_of_lows[j];
      unsigned high_buf[MAX_UP];
      unsigned high_buf_size = 0;
      for (unsigned k = 0; k < verts_per_low; ++k) {
        unsigned vert = up_vert[high_verts_of_low[k]];
        unsigned first_down = lows_of_verts_offsets[vert];
        unsigned end_down = lows_of_verts_offsets[vert + 1];
        if (high_buf_size)
          high_buf_size = intersect(
              high_buf,
              high_buf_size,
              lows_of_verts + first_down,
              end_down - first_down);
        else
          high_buf_size = copy(
              lows_of_verts + first_down,
              high_buf,
              end_down - first_down);
        assert(high_buf_size);
      }
      assert(high_buf_size == 1);
      lows_of_high[j] = high_buf[0];
    }
  }
  return lows_of_highs;
}
