#include "ints.h"
#include <stdlib.h>

void ints_zero(unsigned* a, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    a[i] = 0;
}

void ints_copy(unsigned const* a, unsigned* b, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
}

unsigned ints_max(unsigned const* a, unsigned n)
{
  unsigned max = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

unsigned* ints_exscan(unsigned const* a, unsigned n)
{
  unsigned* o = malloc(sizeof(unsigned) * (n + 1));
  unsigned sum = 0;
  o[0] = 0;
  for (unsigned i = 0; i < n; ++i) {
    sum += a[i];
    o[i + 1] = sum;
  }
  return o;
}

unsigned* ints_unscan(unsigned const* a, unsigned n)
{
  unsigned* o = malloc(sizeof(unsigned) * n);
  for (unsigned i = 0; i < n; ++i)
    o[i] = a[i + 1] - a[i];
  return o;
}

unsigned* ints_negate(unsigned const* a, unsigned n)
{
  unsigned* o = malloc(sizeof(unsigned) * n);
  for (unsigned i = 0; i < n; ++i)
    o[i] = !a[i];
  return o;
}
