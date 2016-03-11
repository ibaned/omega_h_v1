#include "ints.hpp"

#include "loop.hpp"

#if defined(LOOP_CUDA_HPP) || \
      (defined(LOOP_KOKKOS_HPP) && defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA))
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#elif defined(LOOP_OPENMP_HPP) || \
      (defined(LOOP_KOKKOS_HPP) && defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP))
#include <omp.h>
#else
#include <cstdlib>
#endif

#if defined(LOOP_SERIAL_HPP) || \
    defined(LOOP_OPENMP_HPP) || \
    (defined(LOOP_KOKKOS_HPP) && !defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA))
static unsigned* serial_uints_exscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_HOST_MALLOC(unsigned, (n + 1));
  unsigned sum = 0;
  o[0] = 0;
  for (unsigned i = 0; i < n; ++i) {
    sum += a[i];
    o[i + 1] = sum;
  }
  return o;
}
#endif

#if defined(LOOP_CUDA_HPP) || \
      (defined(LOOP_KOKKOS_HPP) && defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA))

unsigned uints_max(unsigned const* a, unsigned n)
{
  unsigned max = 0;
  thrust::device_ptr<unsigned const> p(a);
  max = thrust::reduce(p, p + n, 0, thrust::maximum<unsigned>());
  return max;
}

unsigned* uints_exscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned , n + 1);
  thrust::device_ptr<unsigned const> inp(a);
  thrust::device_ptr<unsigned> outp(o);
  thrust::exclusive_scan(inp, inp + n, outp);
  /* fixup the last element quirk */
  unsigned sum = thrust::reduce(inp, inp + n);
  CUDACALL(cudaMemcpy(o + n, &sum, sizeof(unsigned), cudaMemcpyHostToDevice));
  return o;
}

unsigned uints_sum(unsigned const* a, unsigned n)
{
  unsigned sum = 0;
  thrust::device_ptr< unsigned int const> p (a);
  sum = thrust::reduce(p, p + n, (unsigned)0, thrust::plus<unsigned>());
  return sum;
}

unsigned long ulongs_max(unsigned long const* a, unsigned n)
{
  unsigned long max = 0;
  thrust::device_ptr< unsigned long const> p(a);
  max = thrust::reduce(p, p + n, LONG_MIN, thrust::maximum<unsigned long>());
  return max;
}

#elif defined(LOOP_OPENMP_HPP) || \
      (defined(LOOP_KOKKOS_HPP) && defined(KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP))

unsigned uints_max(unsigned const* a, unsigned n)
{
  unsigned max = 0;
  #pragma omp parallel
  {
    unsigned thread_max = 0;
    #pragma omp for
    for (unsigned i = 0; i < n; ++i)
      if (a[i] > thread_max)
        thread_max = a[i];
    #pragma omp critical
    {
      if (thread_max > max)
        max = thread_max;
    }
  }
  return max;
}

unsigned* uints_exscan(unsigned const* a, unsigned n)
{
  unsigned nthreads = (unsigned) omp_get_max_threads();
  unsigned* thread_sums = LOOP_HOST_MALLOC(unsigned, nthreads);
  #pragma omp parallel
  {
    unsigned thread_sum = 0;
    #pragma omp for schedule(static)
    for (unsigned i = 0; i < n; ++i)
      thread_sum += a[i];
    thread_sums[omp_get_thread_num()] = thread_sum;
  }
  unsigned* thread_exscan = serial_uints_exscan(thread_sums, nthreads);
  loop_host_free(thread_sums);
  unsigned* o = LOOP_HOST_MALLOC(unsigned, (n + 1));
  o[0] = 0;
  #pragma omp parallel
  {
    unsigned sum = thread_exscan[omp_get_thread_num()];
    #pragma omp for schedule(static)
    for (unsigned i = 0; i < n; ++i) {
      sum += a[i];
      o[i + 1] = sum;
    }
  }
  loop_host_free(thread_exscan);
  return o;
}

unsigned uints_sum(unsigned const* a, unsigned n)
{
  unsigned sum = 0;
  #pragma omp parallel for reduction (+:sum)
  for (unsigned i = 0; i < n; ++i)
    sum = sum + a[i];
  return sum;
}

unsigned long ulongs_max(unsigned long const* a, unsigned n)
{
  unsigned long max = 0;
  #pragma omp parallel
  {
    unsigned long thread_max = 0;
    #pragma omp for
    for (unsigned i = 0; i < n; ++i)
      if (a[i] > thread_max)
        thread_max = a[i];
    #pragma omp critical
    {
      if (thread_max > max)
        max = thread_max;
    }
  }
  return max;
}

#else

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
  return serial_uints_exscan(a, n);
}

unsigned uints_sum(unsigned const* a, unsigned n)
{
  unsigned sum = 0;
  for (unsigned i = 0; i < n; ++i)
    sum += a[i];
  return sum;
}

unsigned long ulongs_max(unsigned long const* a, unsigned n)
{
  unsigned long max = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

#endif

#define GENERAL_LINEAR(T, name) \
LOOP_KERNEL(name##_linear_kern, T* out, T stride) \
  out[i] = i * stride; \
} \
T* name##_linear(unsigned n, T stride) \
{ \
  T* out = LOOP_MALLOC(T, n); \
  LOOP_EXEC(name##_linear_kern, n, out, stride); \
  return out; \
}

GENERAL_LINEAR(unsigned, uints)
GENERAL_LINEAR(unsigned long, ulongs)

LOOP_KERNEL(unscan_kern, unsigned const* a, unsigned* o)
  o[i] = a[i + 1] - a[i];
}

unsigned* uints_unscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  LOOP_EXEC(unscan_kern, n, a, o);
  return o;
}

LOOP_KERNEL(negate_kern, unsigned const* a, unsigned* o)
  o[i] = !a[i];
}

unsigned* uints_negate(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  LOOP_EXEC(negate_kern, n, a, o);
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

LOOP_KERNEL(scale_kern, unsigned const* a, unsigned s, unsigned* o)
  o[i] = s * a[i];
}

unsigned* uints_scale(unsigned const* a, unsigned n, unsigned s)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  LOOP_EXEC(scale_kern, n, a, s, o);
  return o;
}
