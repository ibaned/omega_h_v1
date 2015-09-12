#include "star.h"
#include "tables.h"
#include "ints.h"
#include "loop.h"
#include <assert.h>

/* This is the #2 most expensive function, takes up 30% of
   refinement time !
   If you are going to optimize, optimize here !
 */

static unsigned get_ent_star(
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned lows_per_high,
    unsigned low,
    unsigned* star)
{
  unsigned first_up = highs_of_lows_offsets[low];
  unsigned last_up = highs_of_lows_offsets[low + 1];
  unsigned size = 0;
  for (unsigned i = first_up; i < last_up; ++i) {
    unsigned high = highs_of_lows[i];
    for (unsigned j = 0; j < lows_per_high; ++j) {
      unsigned star_low = lows_of_highs[high * lows_per_high + j];
      if (star_low == low)
        continue;
      size = add_unique(star, size, star_low);
    }
  }
  return size;
}

void get_star(
    unsigned low_dim,
    unsigned high_dim,
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned** star_offsets_out,
    unsigned** star_out)
{
  unsigned star_buf[MAX_UP * MAX_DOWN];
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned* degrees = loop_malloc(sizeof(unsigned) * nlows);
  ints_zero(degrees, nlows);
  for (unsigned i = 0; i < nlows; ++i)
    degrees[i] = get_ent_star(
        highs_of_lows_offsets,
        highs_of_lows,
        lows_of_highs,
        lows_per_high,
        i,
        star_buf);
  unsigned* star_offsets = ints_exscan(degrees, nlows);
  loop_free(degrees);
  unsigned sum_degrees = star_offsets[nlows];
  unsigned* star = loop_malloc(sizeof(unsigned) * sum_degrees);
  for (unsigned i = 0; i < nlows; ++i) {
    get_ent_star(
        highs_of_lows_offsets,
        highs_of_lows,
        lows_of_highs,
        lows_per_high,
        i,
        star_buf);
    unsigned first_star = star_offsets[i];
    unsigned last_star = star_offsets[i + 1];
    for (unsigned j = first_star; j < last_star; ++j) {
      assert(j >= first_star);
      star[j] = star_buf[j - first_star];
    }
  }
  *star_offsets_out = star_offsets;
  *star_out = star;
}
