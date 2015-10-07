#include "up_from_down.h"

#include "ints.h"
#include "loop.h"
#include "tables.h"

void up_from_down(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned nlows,
    unsigned const* lows_of_highs,
    unsigned** offsets_out,
    unsigned** highs_of_lows_out,
    unsigned** directions_out)
{
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned nuses = nhighs * lows_per_high;
  unsigned* degrees = loop_malloc(sizeof(unsigned) * nlows);
  uints_zero(degrees, nlows);
  for (unsigned i = 0; i < nuses; ++i)
    degrees[lows_of_highs[i]]++;
  unsigned* offsets = uints_exscan(degrees, nlows);
  unsigned* highs_of_lows = loop_malloc(sizeof(unsigned) * nuses);
  unsigned* directions = 0;
  if (directions_out)
    directions = loop_malloc(sizeof(unsigned) * nuses);
  uints_zero(degrees, nlows);
  for (unsigned i = 0; i < nuses; ++i) {
    unsigned high = i / lows_per_high;
    unsigned direction = i % lows_per_high;
    unsigned down = lows_of_highs[i];
    unsigned o = offsets[down];
    unsigned j = degrees[down]++;
    highs_of_lows[o + j] = high;
    if (directions_out)
      directions[o + j] = direction;
  }
  loop_free(degrees);
  *offsets_out = offsets;
  *highs_of_lows_out = highs_of_lows;
  if (directions_out)
    *directions_out = directions;
}
