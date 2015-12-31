#ifndef LOOP_CUDA_H
#define LOOP_CUDA_H

#include "loop_host.h"

#include <assert.h>

void* loop_cuda_malloc(unsigned long n);
#define LOOP_CUDA_MALLOC(T, n) \
  ((T*)loop_cuda_malloc(sizeof(T) * (n)))
void loop_cuda_free(void* p);

void* loop_cuda_to_host(void const* p, unsigned long n);
void* loop_cuda_to_device(void const* p, unsigned long n);

void loop_cuda_memcpy(void* dst, void const* src, unsigned long n);

static inline __device__ unsigned loop_cuda_atomic_increment(unsigned* p)
{
  return atomicAdd(p, 1);
}

#define loop_to_host loop_cuda_to_host
#define loop_to_device loop_cuda_to_device
#define loop_memcpy loop_cuda_memcpy

#define loop_atomic_increment loop_cuda_atomic_increment

#define LOOP_MALLOC(T, n) LOOP_CUDA_MALLOC(T, n)
#define loop_free loop_cuda_free

#define LOOP_KERNEL(fname, ...) \
static __global__ void fname(unsigned n, __VA_ARGS__) \
{ \
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;\
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
fname<<< loop_ceildiv((n), LOOP_BLOCK_SIZE), LOOP_BLOCK_SIZE >>>(n,__VA_ARGS__);\
CUDACALL(cudaGetLastError()); \
} while(0)

unsigned loop_size(void);

#define LOOP_IN __device__
#define LOOP_INOUT __device__ __host__

#endif
