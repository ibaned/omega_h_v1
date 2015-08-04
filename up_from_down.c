#include "up_from_down.h"
#include "tables.h"
#include "ints.h"
#include <stdlib.h>

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
  unsigned* degrees = malloc(sizeof(unsigned) * nlows);
  ints_zero(degrees, nlows);
  for (unsigned i = 0; i < nuses; ++i)
    degrees[lows_of_highs[i]]++;
  unsigned* offsets = ints_exscan(degrees, nlows);
  unsigned* highs_of_lows = malloc(sizeof(unsigned) * nuses);
  unsigned* directions = 0;
  if (directions_out)
    directions = malloc(sizeof(unsigned) * nuses);
  ints_zero(degrees, nlows);
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
  free(degrees);
  *offsets_out = offsets;
  *highs_of_lows_out = highs_of_lows;
  if (directions_out)
    *directions_out = directions;
}
