#include "subset.h"
#include "ints.h"
#include <stdlib.h>

unsigned* ints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* offsets)
{
  unsigned nsub = offsets[n];
  unsigned* out = malloc(sizeof(unsigned) * nsub * width);
  for (unsigned i = 0; i < n; ++i) {
    if (offsets[i] == offsets[i + 1])
      continue;
    for (unsigned j = 0; j < width; ++j)
      out[offsets[i] * width + j] = a[i * width + j];
  }
  return out;
}
