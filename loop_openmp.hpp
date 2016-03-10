#ifndef LOOP_OPENMP_HPP
#define LOOP_OPENMP_HPP

#include <assert.h>

#include "loop_host.hpp"

#define LOOP_MALLOC(T, n) LOOP_HOST_MALLOC(T, n)
#define loop_free loop_host_free

static inline unsigned loop_openmp_atomic_increment(unsigned* p)
{
  unsigned o;
#pragma omp atomic capture
  o = (*p)++;
  return o;
}

#define loop_atomic_increment loop_openmp_atomic_increment

unsigned loop_size(void);

#define LOOP_IN
#define LOOP_INOUT

#define LOOP_CONST

#define LOOP_NORETURN(x) assert(0)

#endif
