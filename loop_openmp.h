#ifndef LOOP_OPENMP_H
#define LOOP_OPENMP_H

#include "loop_host.h"

#define LOOP_MALLOC(T, n) LOOP_HOST_MALLOC(T, n)
#define LOOP_FREE(p) LOOP_HOST_FREE(p)

#define LOOP_KERNEL(fname, ...) \
static void fname(__VA_ARGS__, unsigned i) \
{

#define LOOP_EXEC(fname, n, ...) \
_Pragma("omp parallel for") \
for (unsigned i = 0; i < n; ++i) \
  fname(__VA_ARGS__, i);

#endif
