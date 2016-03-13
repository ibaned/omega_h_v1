#include "arrays.hpp"

#include <algorithm> //for std::max

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

template <typename T>
T* array_to_device(T const* a, unsigned n)
{
#if defined(LOOP_CUDA_HPP) || \
    (defined(LOOP_KOKKOS_HPP) && defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA))
  T* b = LOOP_MALLOC(T, n);
  CUDACALL(cudaMemcpy(b, a, n * sizeof(T), cudaMemcpyHostToDevice));
  return b;
#else
  return copy_array<T>(a, n);
#endif
}

template unsigned char* array_to_device(unsigned char const* a, unsigned n);
template unsigned* array_to_device(unsigned const* a, unsigned n);
template unsigned long* array_to_device(unsigned long const* a, unsigned n);
template double* array_to_device(double const* a, unsigned n);

template <typename T>
T* array_to_host(T const* a, unsigned n)
{
#if defined(LOOP_CUDA_HPP) || \
    (defined(LOOP_KOKKOS_HPP) && defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA))
  T* b = LOOP_MALLOC(T, n);
  CUDACALL(cudaMemcpy(b, a, n * sizeof(T), cudaMemcpyDeviceToHost));
  return b;
#else
  return copy_array<T>(a, n);
#endif
}

template unsigned char* array_to_host(unsigned char const* a, unsigned n);
template unsigned* array_to_host(unsigned const* a, unsigned n);
template unsigned long* array_to_host(unsigned long const* a, unsigned n);
template double* array_to_host(double const* a, unsigned n);

template <typename T>
LOOP_KERNEL(reorder_kern, T const* a, unsigned width,
    unsigned const* old_to_new, T* o)
  unsigned j = old_to_new[i];
  for (unsigned k = 0; k < width; ++k)
    o[j * width + k] = a[i * width + k];
}
template <typename T>
LOOP_KERNEL(reorder_inv_kern, T const* a, unsigned width,
    unsigned const* new_to_old, T* o)
  unsigned j = new_to_old[i];
  for (unsigned k = 0; k < width; ++k)
    o[i * width + k] = a[j * width + k];
}
template <typename T>
T* reorder_array(T const* a, unsigned const* old_to_new,
    unsigned n, unsigned width)
{
  T* o = LOOP_MALLOC(T, n * width);
  LOOP_EXEC(reorder_kern<T>, n, a, width, old_to_new, o);
  return o;
}
template <typename T>
T* reorder_array_inv(T const* a, unsigned const* new_to_old,
    unsigned n, unsigned width)
{
  T* o = LOOP_MALLOC(T, n * width);
  LOOP_EXEC(reorder_inv_kern<T>, n, a, width, new_to_old, o);
  return o;
}

template unsigned* reorder_array(unsigned const* a,
    unsigned const* old_to_new, unsigned n, unsigned width);
template unsigned long* reorder_array(unsigned long const* a,
    unsigned const* old_to_new, unsigned n, unsigned width);
template double* reorder_array(double const* a,
    unsigned const* old_to_new, unsigned n, unsigned width);

template unsigned* reorder_array_inv(unsigned const* a,
    unsigned const* new_to_old, unsigned n, unsigned width);
template unsigned long* reorder_array_inv(unsigned long const* a,
    unsigned const* new_to_old, unsigned n, unsigned width);
template double* reorder_array_inv(double const* a,
    unsigned const* new_to_old, unsigned n, unsigned width);

#if defined(LOOP_CUDA_HPP) || \
    (defined(LOOP_KOKKOS_HPP) && defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA))
template <typename T>
T array_at(T const* a, unsigned i)
{
  T x;
  CUDACALL(cudaMemcpy(&x, a + i, sizeof(T), cudaMemcpyDeviceToHost));
  return x;
}
#else
template <typename T>
T array_at(T const* a, unsigned i)
{
  return a[i];
}
#endif

template unsigned char array_at(unsigned char const* a, unsigned i);
template unsigned array_at(unsigned const* a, unsigned i);

template <typename T>
LOOP_KERNEL(expand_kern, T const* a, unsigned width,
    unsigned const* offsets, T* out)
  unsigned first = offsets[i];
  unsigned end = offsets[i + 1];
  for (unsigned j = first; j < end; ++j)
    for (unsigned k = 0; k < width; ++k)
      out[j * width + k] = a[i * width + k];
}
template <typename T>
void expand_into(T* out, T const* a,
    unsigned const* offsets, unsigned n, unsigned width)
{
  LOOP_EXEC(expand_kern<T>, n, a, width, offsets, out);
}
template <typename T>
T* expand_array(unsigned n, unsigned width,
    T const* a, unsigned const* offsets)
{
  unsigned nout = array_at(offsets, n);
  T* out = LOOP_MALLOC(T, nout * width);
  expand_into<T>(out, a, offsets, n, width);
  return out;
}

template void expand_into(unsigned* out, unsigned const* a,
    unsigned const* offsets, unsigned n, unsigned width);
template void expand_into(unsigned long* out, unsigned long const* a,
    unsigned const* offsets, unsigned n, unsigned width);

template unsigned char* expand_array(unsigned n, unsigned width,
    unsigned char const* a, unsigned const* offsets);
template unsigned* expand_array(unsigned n, unsigned width,
    unsigned const* a, unsigned const* offsets);
template unsigned long* expand_array(unsigned n, unsigned width,
    unsigned long const* a, unsigned const* offsets);
template double* expand_array(unsigned n, unsigned width,
    double const* a, unsigned const* offsets);

template <typename T>
T* concat_arrays(unsigned width,
    T const* a, unsigned na,
    T const* b, unsigned nb)
{
  T* out = LOOP_MALLOC(T, (na + nb) * width);
  array_memcpy<T>(out, a, na * width);
  array_memcpy<T>(out + (na * width), b, nb * width);
  return out;
}

template unsigned* concat_arrays(unsigned width,
    unsigned const* a, unsigned na,
    unsigned const* b, unsigned nb);
template double* concat_arrays(unsigned width,
    double const* a, unsigned na,
    double const* b, unsigned nb);

template <typename T>
LOOP_KERNEL(fill_kern, T* a, T v)
  a[i] = v;
}
template <typename T>
T* filled_array(unsigned n, T v)
{
  T* a = LOOP_MALLOC(T, n);
  LOOP_EXEC(fill_kern<T>, n, a, v);
  return a;
}

template unsigned char* filled_array(unsigned n, unsigned char v);
template unsigned* filled_array(unsigned n, unsigned v);
template unsigned long* filled_array(unsigned n, unsigned long v);
template double* filled_array(unsigned n, double v);

LOOP_KERNEL(max_doubles_into_kern, double const* a, unsigned width,
    unsigned const* offsets, double* out)
  unsigned first = offsets[i];
  unsigned end = offsets[i + 1];
  if (end == first)
    return;
  for (unsigned k = 0; k < width; ++k)
    out[i * width + k] = a[first * width + k];
  for (unsigned j = first + 1; j < end; ++j)
    for (unsigned k = 0; k < width; ++k)
      out[i * width + k] = std::max(out[i * width + k], a[j * width + k]);
}
void doubles_max_into(unsigned n, unsigned width,
    double const* a, unsigned const* offsets, double* out)
{
  LOOP_EXEC(max_doubles_into_kern, n, a, width, offsets, out);
}
