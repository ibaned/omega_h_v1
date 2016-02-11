#ifndef LOOP_CUDA_H
#define LOOP_CUDA_H

#include "loop_host.h"

#include <assert.h>

#define LOOP_IN __device__
#define LOOP_INOUT __device__ __host__

#define LOOP_CONST __constant__

void* loop_cuda_malloc(unsigned long n);
#define LOOP_CUDA_MALLOC(T, n) \
  ((T*)loop_cuda_malloc(sizeof(T) * (n)))
void loop_cuda_free(void* p);

static inline LOOP_IN unsigned
loop_cuda_atomic_increment(unsigned* p)
{
  return atomicAdd(p, 1);
}

#define loop_atomic_increment loop_cuda_atomic_increment

#define LOOP_MALLOC(T, n) LOOP_CUDA_MALLOC(T, n)
#define loop_free loop_cuda_free

#define LOOP_KERNEL(fname, ...) \
static __global__ void fname(unsigned n, __VA_ARGS__) \
{ \
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x; \
  if (i >= n) return;

#define LOOP_BLOCK_SIZE 256

static inline unsigned loop_ceildiv(unsigned a, unsigned b)
{
  unsigned c = a / b;
  if (a % b)
    ++c;
  return c;
}

#define CUDACALL(f) \
do { \
  cudaError_t ret = (f); \
  assert(ret == cudaSuccess); \
} while (0)

#define LOOP_EXEC(fname, n, ...) \
do { \
  fname<<<loop_ceildiv((n),LOOP_BLOCK_SIZE),LOOP_BLOCK_SIZE>>>(n,__VA_ARGS__); \
  CUDACALL(cudaGetLastError()); \
} while(0)

unsigned loop_size(void);

#define LOOP_NORETURN(x) return x

#endif
