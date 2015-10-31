#include <stdio.h>
#include <stdlib.h>

#include "loop_cuda.h"

void* loop_cuda_malloc(unsigned long n)
{
  void* p;
  CUDACALL(cudaMalloc(&p, n));
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


__device__  unsigned cuda_atomic_increment( unsigned * p )
{
	int a= *p;
	atomicAdd( p , 1);
	return a;
}
