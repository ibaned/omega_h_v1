#ifndef LOOP_OPENMP_H
#define LOOP_OPENMP_H

#include <assert.h>

#include "loop_host.h"

#define LOOP_MALLOC(T, n) LOOP_HOST_MALLOC(T, n)
#define loop_free loop_host_free

/* TODO: these could actually be parallelized across
   threads instead of being left as serial operations */
#define LOOP_TO_HOST LOOP_HOST_COPY

static inline unsigned loop_openmp_atomic_increment(unsigned* p)
{
  unsigned o;
#pragma omp atomic capture
  o = (*p)++;
  return o;
}

#define loop_atomic_increment loop_openmp_atomic_increment

#define LOOP_KERNEL(fname, ...) \
static void fname(__VA_ARGS__, unsigned i) \
{

#define LOOP_EXEC(fname, n, ...) \
_Pragma("omp parallel for") \
for (unsigned i = 0; i < n; ++i) \
  fname(__VA_ARGS__, i);

unsigned loop_size(void);

#define LOOP_IN
#define LOOP_INOUT

#define LOOP_NORETURN(x) assert(0)

#endif
