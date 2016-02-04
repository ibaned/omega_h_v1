#include <stdio.h>
#include <stdlib.h>

#include "loop_cuda.h"

void* loop_cuda_malloc(unsigned long n)
{
  void* p;
  if (!n)
    n = 1;
#if USE_CUDA_MALLOC_MANAGED
  CUDACALL(cudaMallocManaged(&p, n));
#else
  CUDACALL(cudaMalloc(&p, n));
#endif
  return p;
}

void loop_cuda_free(void* p)
{
  CUDACALL(cudaFree(p));
}

void* loop_cuda_to_host(void const* p, unsigned long n)
{
  void* out = loop_host_malloc(n);
  CUDACALL(cudaMemcpy(out, p, n, cudaMemcpyDeviceToHost));
  return out;
}

void* loop_cuda_to_device(void const* p, unsigned long n)
{
  void* out = loop_cuda_malloc(n);
  CUDACALL(cudaMemcpy(out, p, n, cudaMemcpyHostToDevice));
  return out;
}

unsigned loop_size(void)
{
  int device;
  cudaGetDevice(&device);
  /* this may be an over-estimate ? */
  int sm_count;
  cudaDeviceGetAttribute(&sm_count, cudaDevAttrMultiProcessorCount, device);
  int sm_size;
  cudaDeviceGetAttribute(&sm_size, cudaDevAttrMaxThreadsPerMultiProcessor, device);
  return (unsigned) (sm_count * sm_size);
}
