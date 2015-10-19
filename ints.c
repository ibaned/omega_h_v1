#include "ints.h"

#include <stdlib.h>

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

static int uints_less(void const* a, void const* b)
{
  unsigned const* pa = (unsigned const*) a;
  unsigned const* pb = (unsigned const*) b;
  if (pa < pb)
    return -1;
  if (pa > pb)
    return 1;
  return 0;
}

unsigned* uints_sort(unsigned const* a, unsigned n)
{
  unsigned* out = uints_copy(a, n);
  qsort(out, n, sizeof(unsigned), uints_less);
  return out;
}

void uints_unique(unsigned const* a, unsigned n,
    unsigned* nunique, unsigned** unique)
{
  unsigned* sorted = uints_sort(a, n);
  unsigned* jump = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    jump[i] = ((i == 0) || (sorted[i - 1] != sorted[i]));
  unsigned* scan = uints_exscan(jump, n);
  *nunique = scan[n];
  *unique = LOOP_MALLOC(unsigned, *nunique);
  for (unsigned i = 0; i < n; ++i)
    if (jump[i])
      (*unique)[scan[i]] = sorted[i];
  loop_free(sorted);
  loop_free(jump);
  loop_free(scan);
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

unsigned char* uchars_copy(unsigned char const* a, unsigned n)
{
  unsigned char* b = LOOP_MALLOC(unsigned char, n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}
