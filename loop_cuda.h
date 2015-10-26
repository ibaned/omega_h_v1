#ifndef LOOP_CUDA_H
#define LOOP_CUDA_H

#include "loop_host.h"

void* loop_cuda_malloc(unsigned long n);
#define LOOP_CUDA_MALLOC(T, n) \
  ((T*)loop_cuda_malloc(sizeof(T) * (n)))
void loop_cuda_free(void* p);

void* loop_cuda_to_host(void const* p, unsigned long n);
void* loop_cuda_to_device(void const* p, unsigned long n);

#define loop_to_host loop_cuda_to_host
#define loop_to_device loop_cuda_to_device

#define LOOP_MALLOC(T, n) LOOP_CUDA_MALLOC(T, n)
#define loop_free loop_cuda_free

#define LOOP_KERNEL(fname, ...) \
static __global__ void fname(__VA_ARGS__) \
{ \
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;

#define LOOP_BLOCK_SIZE 256

static inline unsigned loop_ceildiv(unsigned a, unsigned b)
{
  unsigned c = a / b;
  if (a % b)
    ++c;
  return c;
}

#define LOOP_EXEC(fname, n, ...) \
fname<<< loop_ceildiv((n), LOOP_BLOCK_SIZE), LOOP_BLOCK_SIZE >>>(__VA_ARGS__);

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

#endif
