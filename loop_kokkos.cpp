#include "loop_kokkos.hpp"

#ifdef KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA

void* loop_kokkos_malloc(unsigned long n)
{
  void* p;
  if (!n)
    n = 1;
  CUDACALL(cudaMalloc(&p, n));
  return p;
}

void loop_kokkos_free(void* p)
{
  CUDACALL(cudaFree(p));
}

#else

void* loop_kokkos_malloc(unsigned long n)
{
  return loop_host_malloc(n);
}

void loop_kokkos_free(void* p)
{
  return loop_host_free(p);
}

#endif
