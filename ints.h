#ifndef INTS_H
#define INTS_H

unsigned uints_max(unsigned const* a, unsigned n);
unsigned* uints_exscan(unsigned const* a, unsigned n);
unsigned* uints_unscan(unsigned const* a, unsigned n);
unsigned* uints_negate(unsigned const* a, unsigned n);
unsigned* uints_negate_offsets(unsigned const* a, unsigned n);
unsigned* uints_filled(unsigned n, unsigned v);
unsigned* uints_linear(unsigned n, unsigned stride);
unsigned uints_sum(unsigned const* a, unsigned n);
unsigned long ulongs_max(unsigned long const* a, unsigned n);

#ifdef __CUDACC__
__device__ __host__
#endif
static inline unsigned has(unsigned const* a, unsigned n, unsigned e)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] == e)
      return 1;
  return 0;
}

#ifdef __CUDACC__
__device__  __host__
#endif
static inline unsigned add_unique(unsigned* a, unsigned n, unsigned e)
{
  if (has(a, n, e))
    return n;
  a[n] = e;
  return n + 1;
}

#endif
