#ifndef LOOP_KOKKOS_HPP
#define LOOP_KOKKOS_HPP

#include "loop_host.hpp"

#ifdef __clang__
#pragma clang system_header
#endif
#include <Kokkos_Core.hpp>

#ifdef KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA
#define LOOP_IN __device__
#define LOOP_CONST __const__
#else
#define LOOP_IN
#define LOOP_CONST
#endif

#define LOOP_INOUT KOKKOS_FUNCTION

void* loop_kokkos_malloc(unsigned long n);
#define LOOP_KOKKOS_MALLOC(T, n) \
  ((T*)loop_kokkos_malloc(sizeof(T) * (n)))
void loop_kokkos_free(void* p);

#define LOOP_MALLOC(T, n) LOOP_KOKKOS_MALLOC(T, n)
#define loop_free loop_kokkos_free

#if defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA )
static inline LOOP_IN unsigned
loop_cuda_atomic_increment(unsigned* p)
{
  return atomicAdd(p, 1);
}
#define loop_atomic_increment loop_cuda_atomic_increment
#elif defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP )
static inline unsigned loop_openmp_atomic_increment(unsigned* p)
{
  unsigned o;
#pragma omp atomic capture
  o = (*p)++;
  return o;
}
#define loop_atomic_increment loop_openmp_atomic_increment
#else
#define loop_atomic_increment loop_host_atomic_increment
#endif

#define LOOP_KERNEL(fname, ...) \
KOKKOS_INLINE_FUNCTION \
static void fname(unsigned i, __VA_ARGS__) \
{

#define LOOP_EXEC(fname, n, ...) \
Kokkos::parallel_for(n, \
    KOKKOS_LAMBDA (unsigned i) { fname(i, __VA_ARGS__); })

#if defined( __CUDACC__ )
#define LOOP_NORETURN(x) return x
#else
#include <cassert>
#define LOOP_NORETURN(x) assert(0)
#endif

#if defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA )
#define CUDACALL(f) \
do { \
  cudaError_t ret = (f); \
  assert(ret == cudaSuccess); \
} while (0)
#endif

#endif
