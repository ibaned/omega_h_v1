#include "up_from_down.hpp"

#include "invert_map.hpp"
#include "loop.hpp"
#include "tables.hpp"

LOOP_KERNEL(separate,
    unsigned* directions,
    unsigned* highs_of_lows,
    unsigned  lows_per_high)

  unsigned both = highs_of_lows[i];
  unsigned high = both/ lows_per_high;
  highs_of_lows[i] = high;
  if(directions)
  {
    unsigned dir = both % lows_per_high;
    directions[i] = dir;
  }
}

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
  unsigned* highs_of_lows;
  unsigned* offsets;
  invert_map(nuses, lows_of_highs, nlows, &highs_of_lows, &offsets);
  unsigned* directions = 0;
  if (directions_out)
    directions = LOOP_MALLOC(unsigned, nuses);
  LOOP_EXEC(separate, nuses, directions, highs_of_lows, lows_per_high);
  *offsets_out = offsets;
  *highs_of_lows_out = highs_of_lows;
  if (directions_out)
    *directions_out = directions;
}
