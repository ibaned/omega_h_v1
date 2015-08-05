#include "concat.h"
#include <stdlib.h>
#include <string.h>

static void* concat_general(
    unsigned typesize,
    unsigned narrs,
    unsigned width,
    unsigned const* sizes,
    void const* const* arrs)
{
  unsigned total_bytes = 0;
  for (unsigned i = 0; i < narrs; ++i)
    total_bytes += sizes[i];
  total_bytes *= width * typesize;
  void* out = malloc(total_bytes);
  void* p = out;
  for (unsigned i = 0; i < narrs; ++i) {
    unsigned arr_bytes = sizes[i] * width * typesize;
    memcpy(p, arrs[i], arr_bytes);
    p += arr_bytes;
  }
  return out;
}

unsigned* concat_ints(
    unsigned narrs,
    unsigned width,
    unsigned const* sizes,
    unsigned const* const* arrs)
{
  return concat_general(sizeof(unsigned), narrs, width, sizes,
      (void const* const*) arrs);
}

double* concat_doubles(
    unsigned narrs,
    unsigned width,
    unsigned const* sizes,
    double const* const* arrs)
{
  return concat_general(sizeof(double), narrs, width, sizes,
      (void const* const*) arrs);
}
