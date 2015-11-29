#include "arrays.h"

#include "loop.h"

#define GENERIC_COPY(T, name) \
LOOP_KERNEL(copy_##name##_kern, T const* a, T* b) \
  b[i] = a[i]; \
} \
T* name##_copy(T const* a, unsigned n) \
{ \
  T* b = LOOP_MALLOC(T, n); \
  LOOP_EXEC(copy_##name##_kern, n, a, b); \
  return b; \
}

GENERIC_COPY(unsigned char, uchars)
GENERIC_COPY(unsigned, uints)
GENERIC_COPY(unsigned long, ulongs)
GENERIC_COPY(double, doubles)

#define GENERIC_SHUFFLE(T, name) \
LOOP_KERNEL(shuffle_##name##_kern, T const* a, unsigned width, \
    unsigned const* out_of_in, T* o) \
  unsigned j = out_of_in[i]; \
  for (unsigned k = 0; k < width; ++k) \
    o[j * width + k] = a[i * width + k]; \
} \
T* name##_shuffle(unsigned n, T const* a, \
    unsigned width, unsigned const* out_of_in) \
{ \
  T* o = LOOP_MALLOC(T, n); \
  LOOP_EXEC(shuffle_##name##_kern, n, a, width, out_of_in, o); \
  return o; \
} \
LOOP_KERNEL(unshuffle_##name##_kern, T const* a, unsigned width, \
    unsigned const* out_of_in, T* o) \
  unsigned j = out_of_in[i]; \
  for (unsigned k = 0; k < width; ++k) \
    o[i * width + k] = a[j * width + k]; \
} \
T* name##_unshuffle(unsigned n, T const* a, \
    unsigned width, unsigned const* out_of_in) \
{ \
  T* o = LOOP_MALLOC(T, n); \
  LOOP_EXEC(unshuffle_##name##_kern, n, a, width, out_of_in, o); \
  return o; \
}

GENERIC_SHUFFLE(unsigned, uints)
GENERIC_SHUFFLE(unsigned long, ulongs)
GENERIC_SHUFFLE(double, doubles)

#define GENERIC_EXPAND(T, name) \
LOOP_KERNEL(expand_##name##_kern, T const* a, unsigned width, \
    unsigned const* offsets, T* o) \
  unsigned first = offsets[i]; \
  unsigned end = offsets[i + 1]; \
  for (unsigned j = first; j < end; ++j) \
    for (unsigned k = 0; k < width; ++k) \
      o[j * width + k] = a[i * width + k]; \
} \
T* name##_expand(unsigned n, T const* a, \
    unsigned width, unsigned const* offsets) \
{ \
  unsigned nout = offsets[n]; \
  T* o = LOOP_MALLOC(T, nout * width); \
  LOOP_EXEC(expand_##name##_kern, n, a, width, offsets, o); \
  return o; \
}

GENERIC_EXPAND(unsigned, uints)
GENERIC_EXPAND(unsigned long, ulongs)
GENERIC_EXPAND(double, doubles)
