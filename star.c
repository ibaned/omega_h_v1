#include "star.h"

#include <assert.h>

#include "ints.h"
#include "loop.h"
#include "tables.h"

/* This is the #2 most expensive function, takes up 30% of
   refinement time !
   If you are going to optimize, optimize here !
 */



#ifdef __CUDACC__
__device__
#endif
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

LOOP_KERNEL( Degree_Shift,
	unsigned *degrees,
	unsigned const* highs_of_lows_offsets,
	unsigned const* highs_of_lows,
	unsigned const* lows_of_highs,
	unsigned lows_per_high,
	unsigned* star_buf)

  degrees[i] = get_ent_star(
	  	highs_of_lows_offsets,
		highs_of_lows,
		lows_of_highs,
		lows_per_high,
		i,
        star_buf);
}

LOOP_KERNEL( Shift,
	unsigned const* highs_of_lows_offsets,
	unsigned const* highs_of_lows,
	unsigned const* lows_of_highs,
	unsigned lows_per_high,
	unsigned* star_buf,
	unsigned* star,
	unsigned* star_offsets)

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
    //assert(j >= first_star);
    star[j] = star_buf[j - first_star];
  }
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
#ifndef NDEBUG
  for (unsigned i = 0; i < MAX_UP * MAX_DOWN; ++i)
    star_buf[i] = INVALID;
#endif
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned* degrees = LOOP_MALLOC(unsigned, nlows);
  uints_zero(degrees, nlows);
  LOOP_EXEC(Degree_Shift, nlows,
    degrees,
    highs_of_lows_offsets,
    highs_of_lows,
    lows_of_highs,
    lows_per_high,
    star_buf);
/*
  for (unsigned i = 0; i < nlows; ++i)
    degrees[i] = get_ent_star(
        highs_of_lows_offsets,
        highs_of_lows,
        lows_of_highs,
        lows_per_high,
        i,
        star_buf);
*/
  unsigned* star_offsets = uints_exscan(degrees, nlows);
  loop_free(degrees);
  unsigned sum_degrees = star_offsets[nlows];
  unsigned* star = LOOP_MALLOC(unsigned, sum_degrees);
  LOOP_EXEC( Shift ,
    nlows ,
	highs_of_lows_offsets,
	highs_of_lows,
	lows_of_highs,
	lows_per_high,
	star_buf,
	star,
	star_offsets);
/*
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
*/
  *star_offsets_out = star_offsets;
  *star_out = star;
}
