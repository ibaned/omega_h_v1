#include "doubles.h"

#include <float.h>

#include "loop.h"

#ifdef __CUDACC__
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>

double doubles_max(double const* a, unsigned n)
{
  double max = 0;
  thrust::device_ptr<double const> p(a);
  /* DBL_MIN is the smallest positive value, for dealing
     with negatives we should actually initialize to -DBL_MAX */
  max = thrust::reduce(p, p + n, -DBL_MAX, thrust::maximum<double>());
  return max;
}

double doubles_min(double const* a, unsigned n)
{
  double min = 0;
  thrust::device_ptr<double const> p(a);
  min = thrust::reduce(p, p + n, DBL_MAX, thrust::minimum<double>());
  return min;
}

double doubles_sum(double const* a, unsigned n)
{
  double sum = 0;
  thrust::device_ptr< double const> p (a);
  sum = thrust::reduce(p, p + n, (double)0, thrust::plus<double>());
  return sum;
}

double* doubles_exscan(double const* a, unsigned n)
{
  double * o = LOOP_MALLOC(double , n + 1);
  thrust::device_ptr<double const> inp(a);
  thrust::device_ptr<double> outp(o);
  thrust::exclusive_scan(inp, inp + n, outp);
  /* fixup the last element quirk */
  double sum = thrust::reduce(inp, inp + n);
  CUDACALL(cudaMemcpy(o + n, &sum, sizeof(double), cudaMemcpyHostToDevice));
  return o;
}

#else

double doubles_max(double const* a, unsigned n)
{
  /* DBL_MIN is the smallest positive value, for dealing
     with negatives we should actually initialize to -DBL_MAX */
  double max = -DBL_MAX;
  for (unsigned i = 1; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

double doubles_min(double const* a, unsigned n)
{
  double min = DBL_MAX;
  for (unsigned i = 1; i < n; ++i)
    if (a[i] < min)
      min = a[i];
  return min;
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
  double* o = LOOP_MALLOC(double, (n + 1));
  double sum = 0;
  o[0] = 0;
  for (unsigned i = 0; i < n; ++i) {
    sum += a[i];
    o[i + 1] = sum;
  }
  return o;
}

#endif

LOOP_KERNEL(axpy_kern,
    double a,
    double const* x,
    double const* y,
    double* out)
  out[i] = a * x[i] + y[i];
}

void doubles_axpy(double a, double const* x, double const* y,
    double* out, unsigned n)
{
  LOOP_EXEC(axpy_kern, n, a, x, y, out);
}
