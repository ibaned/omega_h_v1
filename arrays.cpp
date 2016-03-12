#include "arrays.hpp"

#include "loop.hpp"

#if defined(LOOP_CUDA_HPP) || \
    (defined(LOOP_KOKKOS_HPP) && defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA))
template <typename T>
void array_memcpy(T* dst, T const* src, unsigned n)
{
  CUDACALL(cudaMemcpy(dst, src, n * sizeof(T), cudaMemcpyDeviceToDevice));
}
#else
template <typename T>
LOOP_KERNEL(memcpy_kern, T* dst, T const* src)
  dst[i] = src[i];
}
template <typename T>
void array_memcpy(T* dst, T const* src, unsigned n)
{
  LOOP_EXEC(memcpy_kern<T>, n, dst, src);
}
#endif

template void array_memcpy(unsigned char* dst,
    unsigned char const* src, unsigned n);
template void array_memcpy(unsigned* dst,
    unsigned const* src, unsigned n);
template void array_memcpy(unsigned long* dst,
    unsigned long const* src, unsigned n);
template void array_memcpy(double* dst,
    double const* src, unsigned n);

template <typename T>
T* copy_array(T const* a, unsigned n)
{
  T* b = LOOP_MALLOC(T, n);
  array_memcpy<T>(b, a, n);
  return b;
}

template unsigned char* copy_array(unsigned char const* a, unsigned n);
template unsigned* copy_array(unsigned const* a, unsigned n);
template unsigned long* copy_array(unsigned long const* a, unsigned n);
template double* copy_array(double const* a, unsigned n);

#if defined(LOOP_CUDA_HPP)
#define GENERIC_TO_DEVICE(T, name) \
T* name##_to_device(T const* a, unsigned n) \
{ \
  T* b = LOOP_MALLOC(T, n); \
  CUDACALL(cudaMemcpy(b, a, n * sizeof(T), cudaMemcpyHostToDevice)); \
  return b; \
}
#else
#define GENERIC_TO_DEVICE(T, name) \
T* name##_to_device(T const* a, unsigned n) \
{ \
  return copy_array<T>(a, n); \
}
#endif

GENERIC_TO_DEVICE(unsigned char, uchars)
GENERIC_TO_DEVICE(unsigned, uints)
GENERIC_TO_DEVICE(unsigned long, ulongs)
GENERIC_TO_DEVICE(double, doubles)

#if defined(LOOP_CUDA_HPP)
#define GENERIC_TO_HOST(T, name) \
T* name##_to_host(T const* a, unsigned n) \
{ \
  T* b = LOOP_HOST_MALLOC(T, n); \
  CUDACALL(cudaMemcpy(b, a, n * sizeof(T), cudaMemcpyDeviceToHost)); \
  return b; \
}
#else
#define GENERIC_TO_HOST(T, name) \
T* name##_to_host(T const* a, unsigned n) \
{ \
  return copy_array<T>(a, n); \
}
#endif

GENERIC_TO_HOST(unsigned char, uchars)
GENERIC_TO_HOST(unsigned, uints)
GENERIC_TO_HOST(unsigned long, ulongs)
GENERIC_TO_HOST(double, doubles)

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
  array_memcpy<T>(out, a, na * width); \
  array_memcpy<T>(out + (na * width), b, nb * width); \
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

#ifdef LOOP_CUDA_HPP
#define GENERIC_AT(T, name) \
T name##_at(T const* a, unsigned i) \
{ \
  T x; \
  CUDACALL(cudaMemcpy(&x, a + i, sizeof(T), cudaMemcpyDeviceToHost)); \
  return x; \
}
#else
#define GENERIC_AT(T, name) \
T name##_at(T const* a, unsigned i) \
{ \
  return a[i]; \
}
#endif

GENERIC_AT(unsigned char, uchars)
GENERIC_AT(unsigned, uints)

#define MAX(a, b) ((b) > (a) ? (b) : (a))

#define GENERIC_MAX_INTO(T, name) \
LOOP_KERNEL(max_##name##_into_kern, T const* a, unsigned width, \
    unsigned const* offsets, T* out) \
  unsigned first = offsets[i]; \
  unsigned end = offsets[i + 1]; \
  if (end == first) \
    return; \
  for (unsigned k = 0; k < width; ++k) { \
    out[i * width + k] = a[first * width + k]; \
  } \
  for (unsigned j = first + 1; j < end; ++j) \
    for (unsigned k = 0; k < width; ++k) { \
      out[i * width + k] = MAX(out[i * width + k], a[j * width + k]); \
    } \
} \
void name##_max_into(unsigned n, unsigned width, \
    T const* a, unsigned const* offsets, \
    T* out) \
{ \
  LOOP_EXEC(max_##name##_into_kern, n, a, width, offsets, out); \
}

GENERIC_MAX_INTO(double, doubles)
