#ifndef LOOP_CUDA_H
#define LOOP_CUDA_H

#include "loop_host.h"

void* loop_cuda_malloc(unsigned long n);
#define LOOP_CUDA_MALLOC(T, n) \
  ((T*)loop_cuda_malloc(sizeof(T) * (n)))
void loop_cuda_free(void* p);

#define LOOP_MALLOC(T, n) LOOP_CUDA_MALLOC(T, n)
#define loop_free loop_cuda_free

#define LOOP_KERNEL(fname, ...) \
static __global__ void fname(__VA_ARGS__) \
{ \
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;

#define LOOP_BLOCK_SIZE 256

static inline int loop_ceildiv(int a, int b)
{
  int c = a / b;
  if (a % b)
    ++c;
  return c;
}

#define LOOP_EXEC(fname, n, ...) \
fname<<< loop_ceildiv((n), LOOP_BLOCK_SIZE), LOOP_BLOCK_SIZE >>>(__VA_ARGS__);

#endif

