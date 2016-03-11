#ifndef LOOP_OPENMP_HPP
#define LOOP_OPENMP_HPP

#include <assert.h>

#include "loop_host.hpp"

#define LOOP_MALLOC(T, n) \
  static_cast<T*>(loop_host_malloc(sizeof(T) * (n)))
#define loop_free loop_host_free

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

#define LOOP_CONST

#define LOOP_NORETURN(x) assert(0)

#endif
