#include "reflect_down.h"

#include <assert.h>

#include "ints.h"
#include "loop.h"
#include "tables.h"

static unsigned copy(
    unsigned const* a,
    unsigned* b,
    unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return n;
}

static unsigned copy_except(
    unsigned const* a,
    unsigned* b,
    unsigned n,
    unsigned exclude_this)
{
  unsigned j = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] != exclude_this)
      b[j++] = a[i];
  return j;
}

static unsigned intersect(
    unsigned* a,
    unsigned na,
    unsigned const* b,
    unsigned nb)
{
  unsigned j = 0;
  for (unsigned i = 0; i < na; ++i)
    if (has(b, nb, a[i]))
      a[j++] = a[i];
  return j;
}

/* This is the #1 most expensive function, takes up 50% of
   refinement time !
   If you are going to optimize, optimize here !
 */
static unsigned* reflect_down_general(
    unsigned dual_mode,
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* verts_of_highs,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts)
{
  unsigned verts_per_high = the_down_degrees[high_dim][0];
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned verts_per_low = the_down_degrees[low_dim][0];
  unsigned* lows_of_highs = LOOP_MALLOC(unsigned, nhighs * lows_per_high);
  unsigned const* const* high_verts_of_lows =
    the_canonical_orders[high_dim][low_dim][0];
  for (unsigned i = 0; i < nhighs; ++i) {
    unsigned const* verts_of_high = verts_of_highs + i * verts_per_high;
    unsigned* lows_of_high = lows_of_highs + i * lows_per_high;
    for (unsigned j = 0; j < lows_per_high; ++j) {
      unsigned const* high_verts_of_low = high_verts_of_lows[j];
      unsigned high_buf[MAX_UP];
      unsigned high_buf_size = 0;
      for (unsigned k = 0; k < verts_per_low; ++k) {
        unsigned vert = verts_of_high[high_verts_of_low[k]];
        unsigned first_use = lows_of_verts_offsets[vert];
        unsigned end_use = lows_of_verts_offsets[vert + 1];
        if (k) {
          high_buf_size = intersect(
              high_buf,
              high_buf_size,
              lows_of_verts + first_use,
              end_use - first_use);
        } else if (dual_mode) {
          assert(end_use - first_use <= MAX_UP);
          high_buf_size = copy_except(
              lows_of_verts + first_use,
              high_buf,
              end_use - first_use,
              i);
        } else {
          assert(end_use - first_use <= MAX_UP);
          high_buf_size = copy(
              lows_of_verts + first_use,
              high_buf,
              end_use - first_use);
        }
      }
      assert(high_buf_size <= 1);
      lows_of_high[j] = ((high_buf_size) ? high_buf[0] : INVALID);
    }
  }
  return lows_of_highs;
}

unsigned* reflect_down(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* verts_of_highs,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts)
{
  return reflect_down_general(0, high_dim, low_dim, nhighs, verts_of_highs,
      lows_of_verts_offsets, lows_of_verts);
}

unsigned* get_dual(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts)
{
  return reflect_down_general(1, elem_dim, elem_dim - 1, nelems, verts_of_elems,
      elems_of_verts_offsets, elems_of_verts);
}
