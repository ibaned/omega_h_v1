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
  char* p = out;
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

void* general_subset(
    unsigned typesize,
    unsigned n,
    unsigned width,
    void const* a,
    unsigned const* offsets)
{
  unsigned nsub = offsets[n];
  void* out = malloc(typesize * nsub * width);
  unsigned stride = typesize * width;
  char const* ip = a;
  char* op = out;
  for (unsigned i = 0; i < n; ++i, ip += stride) {
    if (offsets[i] == offsets[i + 1])
      continue;
    memcpy(op, ip, stride);
    op += stride;
  }
  return out;
}

unsigned* ints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* offsets)
{
  return general_subset(sizeof(unsigned), n, width, a, offsets);
}

double* doubles_subset(
    unsigned n,
    unsigned width,
    double const* a,
    unsigned const* offsets)
{
  return general_subset(sizeof(double), n, width, a, offsets);
}
