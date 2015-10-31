#include "loop_host.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

void* loop_host_malloc(unsigned long n)
{
  if (!n)
    return NULL;
  void* p = malloc(n);
  assert(p);
  return p;
}

unsigned loop_host_atomic_increment(unsigned* p)
{
  unsigned a = *p;
  *p = *p + 1;
  return a;
}

void* loop_host_realloc(void* p, unsigned long n)
{
  if (!n) {
    free(p);
    return NULL;
  }
  if (!p)
    return loop_host_malloc(n);
  void* q = realloc(p, n);
  assert(q);
  return q;
}

void loop_host_free(void* p)
{
  free(p);
}

void* loop_host_copy(void const* p, unsigned long n)
{
  void* out = loop_host_malloc(n);
  memcpy(out, p, n);
  return out;
}
