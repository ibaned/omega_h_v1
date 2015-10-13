#include "ints.h"

#include "loop.h"

void uints_zero(unsigned* a, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    a[i] = 0;
}

unsigned* uints_copy(unsigned const* a, unsigned n)
{
  unsigned* b = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}

unsigned uints_max(unsigned const* a, unsigned n)
{
  unsigned max = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

unsigned* uints_exscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, (n + 1));
  unsigned sum = 0;
  o[0] = 0;
  for (unsigned i = 0; i < n; ++i) {
    sum += a[i];
    o[i + 1] = sum;
  }
  return o;
}

unsigned* uints_unscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    o[i] = a[i + 1] - a[i];
  return o;
}

unsigned* uints_negate(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    o[i] = !a[i];
  return o;
}

unsigned* uints_negate_offsets(unsigned const* a, unsigned n)
{
  unsigned* unscanned = uints_unscan(a, n);
  unsigned* negated = uints_negate(unscanned, n);
  loop_free(unscanned);
  unsigned* out = uints_exscan(negated, n);
  loop_free(negated);
  return out;
}

void uints_fill(unsigned* a, unsigned n, unsigned v)
{
  for (unsigned i = 0; i < n; ++i)
    a[i] = v;
}

unsigned uints_sum(unsigned const* a, unsigned n)
{
  unsigned sum = 0;
  for (unsigned i = 0; i < n; ++i)
    sum += a[i];
  return sum;
}

unsigned long* ulongs_copy(unsigned long const* a, unsigned n)
{
  unsigned long* b = LOOP_MALLOC(unsigned long, n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}

unsigned long ulongs_max(unsigned long const* a, unsigned n)
{
  unsigned long max = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}
