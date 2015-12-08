#include "loop_host.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#if MEASURE_MEMORY
static unsigned long nbytes_allocted;

static void* measure_malloc(unsigned long n)
{
  unsigned long* p = malloc(n + sizeof(unsigned long));
  p[0] = n;
  nbytes_allocted += n;
  return p + 1;
}

static void* measure_realloc(void* p, unsigned long n)
{
  unsigned long* a = p;
  unsigned long on = a[-1];
  void* q = measure_malloc(n);
  memcpy(q, p, on);
  measure_free(p);
  return q;
}

static void measure_free(void* p)
{
  unsigned long* a = p;
  nbytes_allocted -= a[-1];
  free(a - 1);
}
#else
static void* measure_malloc(unsigned long n)
{
  return malloc(n);
}

static void* measure_realloc(void* p, unsigned long n)
{
  return realloc(p, n);
}

static void measure_free(void* p)
{
  free(p);
}
#endif

void* loop_host_malloc(unsigned long n)
{
  if (!n)
    n = 1;
  void* p = measure_malloc(n);
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
  if (!n)
    n = 1;
  void* q = measure_realloc(p, n);
  assert(q);
  return q;
}

void loop_host_free(void* p)
{
  measure_free(p);
}

void* loop_host_copy(void const* p, unsigned long n)
{
  void* out = loop_host_malloc(n);
  memcpy(out, p, n);
  return out;
}
