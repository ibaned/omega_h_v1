#include "ints.h"

#include "loop.h"

#ifdef __CUDACC__
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#else
#include <stdlib.h>
#endif

#ifdef __CUDACC__

unsigned uints_max(unsigned const* a, unsigned n)
{
  unsigned max = 0;
  thrust::device_ptr<unsigned const> p(a);
  max = thrust::reduce(p, p + n, 0, thrust::maximum<unsigned>());
  return max;
}

unsigned* uints_exscan(unsigned const* a, unsigned n)
{
  unsigned * o = LOOP_MALLOC(unsigned , n + 1);
  thrust::device_ptr<unsigned const> inp(a);
  thrust::device_ptr<unsigned> outp(o);
  thrust::exclusive_scan(inp, inp + n, outp);
  /* fixup the last element quirk */
  unsigned sum = thrust::reduce(inp, inp + n);
  CUDACALL(cudaMemcpy(o + n, &sum, sizeof(unsigned), cudaMemcpyHostToDevice));
  return o;
}

unsigned uints_sum(unsigned const* a, unsigned n)
{
  unsigned sum = 0;
  thrust::device_ptr< unsigned int const> p (a);
  sum = thrust::reduce(p, p + n, (unsigned)0, thrust::plus<unsigned>());
  return sum;
}

unsigned long ulongs_max(unsigned long const* a, unsigned n)
{
  unsigned long max = 0;
  thrust::device_ptr< unsigned long const> p(a);
  max = thrust::reduce(p, p + n, LONG_MIN, thrust::maximum<unsigned long>());
  return max;
}

#else

unsigned uints_max(unsigned const* a, unsigned n)
{
  unsigned max = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

unsigned* uints_exscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, (n + 1));
  unsigned sum = 0;
  o[0] = 0;
  for (unsigned i = 0; i < n; ++i) {
    sum += a[i];
    o[i + 1] = sum;
  }
  return o;
}

unsigned uints_sum(unsigned const* a, unsigned n)
{
  unsigned sum = 0;
  for (unsigned i = 0; i < n; ++i)
    sum += a[i];
  return sum;
}

unsigned long ulongs_max(unsigned long const* a, unsigned n)
{
  unsigned long max = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

#endif

unsigned* uints_linear(unsigned n, unsigned stride)
{
  unsigned* scrap = uints_filled(n, stride);
  unsigned* linear = uints_exscan(scrap, n);
  loop_free(scrap);
  return linear;
}

LOOP_KERNEL(unscan_kern, unsigned const* a, unsigned* o)
  o[i] = a[i + 1] - a[i];
}

unsigned* uints_unscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  LOOP_EXEC(unscan_kern, n, a, o);
  return o;
}

LOOP_KERNEL(negate_kern, unsigned const* a, unsigned* o)
  o[i] = !a[i];
}

unsigned* uints_negate(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  LOOP_EXEC(negate_kern, n, a, o);
  return o;
}

unsigned* uints_negate_offsets(unsigned const* a, unsigned n)
{
  unsigned* unscanned = uints_unscan(a, n);
  unsigned* negated = uints_negate(unscanned, n);
  loop_free(unscanned);
  unsigned* out = uints_exscan(negated, n);
  loop_free(negated);
  return out;
}

LOOP_KERNEL(fill_kern, unsigned* a, unsigned v)
  a[i] = v;
}

unsigned* uints_filled(unsigned n, unsigned v)
{
  unsigned* a = LOOP_MALLOC(unsigned, n);
  LOOP_EXEC(fill_kern, n, a, v);
  return a;
}

LOOP_KERNEL(scale_kern, unsigned const* a, unsigned s, unsigned* o)
  o[i] = s * a[i];
}

unsigned* uints_scale(unsigned const* a, unsigned n, unsigned s)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  LOOP_EXEC(scale_kern, n, a, s, o);
  return o;
}
