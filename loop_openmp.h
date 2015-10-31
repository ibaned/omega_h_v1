#ifndef LOOP_OPENMP_H
#define LOOP_OPENMP_H

#include "loop_host.h"

#define LOOP_MALLOC(T, n) LOOP_HOST_MALLOC(T, n)
#define loop_free loop_host_free

#define loop_to_host loop_host_copy
#define loop_to_device loop_host_copy

static inline unsigned loop_openmp_atomic_increment(unsigned* p)
{
  unsigned a = *p;
#pragma omp atomic update
    *p = *p +1;
  return a;
}

#define loop_atomic_increment loop_openmp_atomic_increment

#define LOOP_KERNEL(fname, ...) \
static void fname(__VA_ARGS__, unsigned i) \
{

#define LOOP_EXEC(fname, n, ...) \
_Pragma("omp parallel for") \
for (unsigned i = 0; i < n; ++i) \
  fname(__VA_ARGS__, i);

#endif
