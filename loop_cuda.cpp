#include "loop_cuda.hpp"

void* loop_cuda_malloc(unsigned long n)
{
  void* p;
  if (!n)
    n = 1;
  CUDACALL(cudaMalloc(&p, n));
  return p;
}

void loop_cuda_free(void* p)
{
  CUDACALL(cudaFree(p));
}
