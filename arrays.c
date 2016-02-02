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
  T* o = LOOP_MALLOC(T, n * width); \
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
  T* o = LOOP_MALLOC(T, n * width); \
  LOOP_EXEC(unshuffle_##name##_kern, n, a, width, out_of_in, o); \
  return o; \
}

GENERIC_SHUFFLE(unsigned, uints)
GENERIC_SHUFFLE(unsigned long, ulongs)
GENERIC_SHUFFLE(double, doubles)

#define GENERIC_EXPAND(T, name) \
LOOP_KERNEL(expand_##name##_kern, T const* a, unsigned width, \
    unsigned const* offsets, T* out) \
  unsigned first = offsets[i]; \
  unsigned end = offsets[i + 1]; \
  for (unsigned j = first; j < end; ++j) \
    for (unsigned k = 0; k < width; ++k) \
      out[j * width + k] = a[i * width + k]; \
} \
void name##_expand_into(unsigned n, unsigned width, \
    T const* a, unsigned const* offsets, \
    T* out) \
{ \
  LOOP_EXEC(expand_##name##_kern, n, a, width, offsets, out); \
} \
T* name##_expand(unsigned n, unsigned width, \
    T const* a, unsigned const* offsets) \
{ \
  unsigned nout = uints_at(offsets, n); \
  T* out = LOOP_MALLOC(T, nout * width); \
  name##_expand_into(n, width, a, offsets, out); \
  return out; \
}

GENERIC_EXPAND(unsigned char, uchars)
GENERIC_EXPAND(unsigned, uints)
GENERIC_EXPAND(unsigned long, ulongs)
GENERIC_EXPAND(double, doubles)

#define GENERIC_CONCAT(T, name) \
T* concat_##name(unsigned width, \
    T const* a, unsigned na, \
    T const* b, unsigned nb) \
{ \
  T* out = LOOP_MALLOC(T, (na + nb) * width); \
  LOOP_MEMCPY(T, out, a, na * width); \
  LOOP_MEMCPY(T, out + (na * width), b, nb * width); \
  return out; \
}

GENERIC_CONCAT(unsigned, uints)
GENERIC_CONCAT(double, doubles)

#define GENERIC_FILL(T, name) \
LOOP_KERNEL(name##_fill_kern, T* a, T v) \
  a[i] = v; \
} \
T* name##_filled(unsigned n, T v) \
{ \
  T* a = LOOP_MALLOC(T, n); \
  LOOP_EXEC(name##_fill_kern, n, a, v); \
  return a; \
}

GENERIC_FILL(unsigned char, uchars)
GENERIC_FILL(unsigned, uints)
GENERIC_FILL(unsigned long, ulongs)
GENERIC_FILL(double, doubles)

#ifdef LOOP_CUDA_H
#define GENERIC_AT(T, name) \
T uints_at(T const* a, unsigned i) \
{ \
  T x; \
  CUDACALL(cudaMemcpy(&x, a + i, sizeof(T), cudaMemcpyDeviceToHost)); \
  return x; \
}
#else
#define GENERIC_AT(T, name) \
T uints_at(T const* a, unsigned i) \
{ \
  return a[i]; \
}
#endif

GENERIC_AT(unsigned, uints)
