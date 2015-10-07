#include "doubles.h"

#include "loop.h"

double doubles_max(double const* a, unsigned n)
{
  double max = a[0];
  for (unsigned i = 1; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

double doubles_min(double const* a, unsigned n)
{
  double min = a[0];
  for (unsigned i = 1; i < n; ++i)
    if (a[i] < min)
      min = a[i];
  return min;
}

void doubles_axpy(double a, double const* x, double const* y,
    double* out, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    out[i] = a * x[i] + y[i];
}

double* doubles_copy(double const* a, unsigned n)
{
  double* b = loop_malloc(sizeof(double) * n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}

/* ambitious note to self: this could be one source
   of partitionig/ordering dependence from inputs
   to outputs. */
double doubles_sum(double const* a, unsigned n)
{
  double s = 0;
  for (unsigned i = 0; i < n; ++i)
    s += a[i];
  return s;
}

double* doubles_exscan(double const* a, unsigned n)
{
  double* o = loop_malloc(sizeof(double) * (n + 1));
  double sum = 0;
  o[0] = 0;
  for (unsigned i = 0; i < n; ++i) {
    sum += a[i];
    o[i + 1] = sum;
  }
  return o;
}
