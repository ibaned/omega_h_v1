#ifndef INTS_H
#define INTS_H

unsigned* uints_copy(unsigned const* a, unsigned n);
unsigned uints_max(unsigned const* a, unsigned n);
unsigned* uints_exscan(unsigned const* a, unsigned n);
unsigned* uints_unscan(unsigned const* a, unsigned n);
unsigned* uints_negate(unsigned const* a, unsigned n);
unsigned* uints_negate_offsets(unsigned const* a, unsigned n);
unsigned* uints_filled(unsigned n, unsigned v);
unsigned* uints_linear(unsigned n);
unsigned uints_sum(unsigned const* a, unsigned n);
unsigned* uints_expand(unsigned n, unsigned const* a,
    unsigned width, unsigned const* offsets);

static inline unsigned has(unsigned const* a, unsigned n, unsigned e)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] == e)
      return 1;
  return 0;
}

static inline unsigned add_unique(unsigned* a, unsigned n, unsigned e)
{
  if (has(a, n, e))
    return n;
  a[n] = e;
  return n + 1;
}

unsigned long* ulongs_copy(unsigned long const* a, unsigned n);
unsigned long ulongs_max(unsigned long const* a, unsigned n);

unsigned char* uchars_copy(unsigned char const* a, unsigned n);

#endif
