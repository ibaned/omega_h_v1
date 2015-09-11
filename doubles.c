#include "doubles.h"
#include <stdlib.h>

double doubles_max(double const a[], unsigned n)
{
  double max = a[0];
  for (unsigned i = 1; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

double doubles_min(double const a[], unsigned n)
{
  double min = a[0];
  for (unsigned i = 1; i < n; ++i)
    if (a[i] < min)
      min = a[i];
  return min;
}

void doubles_axpy(double a, double const x[], double const y[],
    double out[], unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    out[i] = a * x[i] + y[i];
}

double* doubles_copy(double const a[], unsigned n)
{
  double* b = malloc(sizeof(double) * n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}
