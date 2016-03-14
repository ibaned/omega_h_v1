#ifndef INTS_HPP
#define INTS_HPP

#include "loop.hpp"

namespace omega_h {

unsigned uints_max(unsigned const* a, unsigned n);
unsigned* uints_exscan(unsigned const* a, unsigned n);
unsigned* uints_unscan(unsigned const* a, unsigned n);
unsigned* uints_negate(unsigned const* a, unsigned n);
unsigned* uints_negate_offsets(unsigned const* a, unsigned n);
unsigned* uints_linear(unsigned n, unsigned stride);
unsigned long* ulongs_linear(unsigned n, unsigned long stride);
unsigned uints_sum(unsigned const* a, unsigned n);
unsigned long ulongs_max(unsigned long const* a, unsigned n);
unsigned* uints_scale(unsigned const* a, unsigned n, unsigned s);

LOOP_INOUT static inline unsigned
has(unsigned const* a, unsigned n, unsigned e)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] == e)
      return 1;
  return 0;
}

LOOP_INOUT static inline unsigned
add_unique(unsigned* a, unsigned n, unsigned e)
{
  if (has(a, n, e))
    return n;
  a[n] = e;
  return n + 1;
}

}

#endif
