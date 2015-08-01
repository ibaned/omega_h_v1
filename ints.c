#include "ints.h"
#include <stdlib.h>

void ints_zero(unsigned* a, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    a[i] = 0;
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
