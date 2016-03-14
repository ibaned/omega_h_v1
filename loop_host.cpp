#include "loop_host.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>

namespace omega_h {

#if MEASURE_MEMORY
static unsigned long memory_usage = 0;
static unsigned long high_water = 0;

static void* measure_malloc(unsigned long n)
{
  unsigned long* p = malloc(n + sizeof(unsigned long));
  p[0] = n;
  memory_usage += n;
  if (memory_usage > high_water)
    high_water = memory_usage;
  return &p[1];
}

static void measure_free(void* p)
{
  if (!p)
    return;
  unsigned long* a = p;
  memory_usage -= a[-1];
  free(&a[-1]);
}

static void* measure_realloc(void* p, unsigned long n)
{
  void* q = measure_malloc(n);
  if (p) {
    unsigned long* a = p;
    unsigned long common = a[-1];
    if (n < common)
      common = n;
    memcpy(q, p, common);
    measure_free(p);
  }
  return q;
}

unsigned long loop_host_memory(void)
{
  return memory_usage;
}

unsigned long loop_host_high_water(void)
{
  return high_water;
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

unsigned long loop_host_memory(void)
{
  return 0;
}

unsigned long loop_host_high_water(void)
{
  return 0;
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

}
