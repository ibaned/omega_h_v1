#include "reflect_down.h"
#include "intersect.h"
#include "tables.h"
#include <stdlib.h>
#include <assert.h>

static unsigned copy_except(
    unsigned const a[],
    unsigned b[],
    unsigned n,
    unsigned exclude_this)
{
  unsigned j = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] != exclude_this)
      b[j++] = a[i];
  return j;
}

static unsigned has(
    unsigned const a[],
    unsigned n,
    unsigned e)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] == e)
      return 1;
  return 0;
}

static unsigned intersect(
    unsigned a[],
    unsigned na,
    unsigned const b[],
    unsigned nb)
{
  unsigned j = 0;
  for (unsigned i = 0; i < na; ++i)
    if (has(b, nb, a[i]))
      a[j++] = a[i];
  return j;
}

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
    unsigned const* verts_of_high = verts_of_highs + i * verts_per_high;
    unsigned* lows_of_high = lows_of_highs + i * lows_per_high;
    unsigned exclude_high = ((high_dim == low_dim) ? i : INVALID);
    for (unsigned j = 0; j < lows_per_high; ++j) {
      unsigned const* high_verts_of_low = high_verts_of_lows[j];
      unsigned high_buf[MAX_UP];
      unsigned high_buf_size = 0;
      for (unsigned k = 0; k < verts_per_low; ++k) {
        unsigned vert = verts_of_high[high_verts_of_low[k]];
        unsigned first_use = lows_of_verts_offsets[vert];
        unsigned end_use = lows_of_verts_offsets[vert + 1];
        if (high_buf_size)
          high_buf_size = intersect(
              high_buf,
              high_buf_size,
              lows_of_verts + first_use,
              end_use - first_use);
        else
          high_buf_size = copy_except(
              lows_of_verts + first_use,
              high_buf,
              end_use - first_use,
              exclude_high);
        assert(high_buf_size);
      }
      assert(high_buf_size == 1);
      lows_of_high[j] = high_buf[0];
    }
  }
  return lows_of_highs;
}
