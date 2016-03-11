#include "loop_cuda.hpp"

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
