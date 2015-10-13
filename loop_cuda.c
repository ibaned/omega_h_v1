#include <stdio.h>
#include <stdlib.h>

#include "loop_cuda.h"

#define CUDACALL(f) \
do { \
  cudaError_t err = f; \
  if (err != cudaSuccess) { \
    const char* errs = cudaGetErrorString(err); \
    fprintf(stderr, "call %s failed at %s:%d : %s\n", \
                    #f, __FILE__, __LINE__, errs); \
    abort(); \
  } \
} while (0)

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
