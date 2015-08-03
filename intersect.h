#ifndef INTERSECT_H
#define INTERSECT_H

static unsigned copy(
    unsigned const a[],
    unsigned b[],
    unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return n;
}

static unsigned has(
    unsigned const a[],
    unsigned n,
    unsigned e)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] == e)
      return 1;
  return 0;
}

static unsigned intersect(
    unsigned a[],
    unsigned na,
    unsigned const b[],
    unsigned nb)
{
  unsigned j = 0;
  for (unsigned i = 0; i < na; ++i)
    if (has(b, nb, a[i]))
      a[j++] = a[i];
  return j;
}

#endif
