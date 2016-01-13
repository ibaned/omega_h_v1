#ifndef LOOP_CUDA_H
#define LOOP_CUDA_H

#include "loop_host.h"

#include <assert.h>

#define LOOP_IN __device__
#define LOOP_INOUT __device__ __host__

void* loop_cuda_malloc(unsigned long n);
#define LOOP_CUDA_MALLOC(T, n) \
  ((T*)loop_cuda_malloc(sizeof(T) * (n)))
void loop_cuda_free(void* p);

void* loop_cuda_to_host(void const* p, unsigned long n);
#define LOOP_CUDA_TO_HOST(T, p, n) \
  ((T*)loop_cuda_to_host(p, sizeof(T) * (n)))
void* loop_cuda_to_device(void const* p, unsigned long n);
#define LOOP_CUDA_TO_DEVICE(T, p, n) \
  ((T*)loop_cuda_to_device(p, sizeof(T) * (n)))
void loop_cuda_memcpy(void* dst, void const* src, unsigned long n);
#define LOOP_CUDA_MEMCPY(T, dst, src, n) \
  loop_cuda_memcpy(dst, src, sizeof(T) * (n))

static inline LOOP_IN unsigned
loop_cuda_atomic_increment(unsigned* p)
{
  return atomicAdd(p, 1);
}

#define LOOP_TO_HOST LOOP_CUDA_TO_HOST
#define LOOP_TO_DEVICE LOOP_CUDA_TO_DEVICE
#define LOOP_MEMCPY LOOP_CUDA_MEMCPY

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
  fname<<<loop_ceildiv((n),LOOP_BLOCK_SIZE),LOOP_BLOCK_SIZE>>>(n,__VA_ARGS__); \
  CUDACALL(cudaGetLastError()); \
} while(0)

unsigned loop_size(void);

#endif
