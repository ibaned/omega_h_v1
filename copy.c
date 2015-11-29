#include "ints.h"
#include "doubles.h"

/* TODO: the split between ints.c and
   doubles.c is not ideal for the most
   basic functions like this copy. */

#include "loop.h"

#define GENERIC_COPY(T, name) \
LOOP_KERNEL(copy_##name##_kern, T const* a, T* b) \
  b[i] = a[i]; \
} \
T* name##_copy(T const* a, unsigned n) \
{ \
  T* b = LOOP_MALLOC(T, n); \
  LOOP_EXEC(copy_##name##_kern, n, a, b); \
  return b; \
}

GENERIC_COPY(unsigned char, uchars)
GENERIC_COPY(unsigned, uints)
GENERIC_COPY(unsigned long, ulongs)
GENERIC_COPY(double, doubles)
