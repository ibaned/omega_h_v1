#include "star.h"

#include <assert.h>
#include <stdio.h>

#include "arrays.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

static LOOP_IN unsigned get_ent_star_general(
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

LOOP_KERNEL(count_general,
    unsigned *degrees,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned lows_per_high)

  unsigned star_buf[MAX_UP * MAX_DOWN];
  degrees[i] = get_ent_star_general(
      highs_of_lows_offsets,
      highs_of_lows,
      lows_of_highs,
      lows_per_high,
      i,
      star_buf);
}

LOOP_KERNEL(fill_general,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned lows_per_high,
    unsigned* star,
    unsigned const* star_offsets)

  get_ent_star_general(
      highs_of_lows_offsets,
      highs_of_lows,
      lows_of_highs,
      lows_per_high,
      i,
      star + star_offsets[i]);
}

static void get_star_general(
    unsigned low_dim,
    unsigned high_dim,
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned** star_offsets_out,
    unsigned** star_out)
{
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned* degrees = uints_filled(nlows, 0);
  LOOP_EXEC(count_general,
    nlows,
    degrees,
    highs_of_lows_offsets,
    highs_of_lows,
    lows_of_highs,
    lows_per_high);
  unsigned* star_offsets = uints_exscan(degrees, nlows);
  loop_free(degrees);
  unsigned sum_degrees = star_offsets[nlows];
  unsigned* star = LOOP_MALLOC(unsigned, sum_degrees);
  LOOP_EXEC(fill_general,
    nlows,
    highs_of_lows_offsets,
    highs_of_lows,
    lows_of_highs,
    lows_per_high,
    star,
    star_offsets);
  *star_offsets_out = star_offsets;
  *star_out = star;
}

void mesh_get_star_general(
    struct mesh* m,
    unsigned low_dim,
    unsigned high_dim,
    unsigned** p_star_offsets,
    unsigned** p_star)
{
  get_star_general(low_dim, high_dim, mesh_count(m, low_dim),
      mesh_ask_up(m, low_dim, high_dim)->offsets,
      mesh_ask_up(m, low_dim, high_dim)->adj,
      mesh_ask_down(m, high_dim, low_dim),
      p_star_offsets, p_star);
}

void mesh_get_star(
    struct mesh* m,
    unsigned low_dim,
    unsigned high_dim,
    unsigned** p_star_offsets,
    unsigned** p_star)
{
  if (low_dim == 0 && high_dim == mesh_dim(m) && mesh_has_dim(m, 1))
    mesh_get_star(m, low_dim, 1, p_star_offsets, p_star);
  else
    mesh_get_star_general(m, low_dim, high_dim, p_star_offsets, p_star);
}
