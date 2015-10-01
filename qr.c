#include "qr.h"
#include <math.h>

static double sign(double x)
{
  return (x < 0) ? -1 : 1;
}

static double square(double x)
{
  return x * x;
}

static unsigned get_reflector(double a[3][3], double v[3], unsigned k)
{
  double cnorm = 0;
  for (unsigned i = k; i < 3; ++i)
    cnorm += square(a[i][k]);
  cnorm = sqrt(cnorm);
  if (cnorm < 1e-14)
    return 0;
  for (unsigned i = 0; i < k; ++i)
    v[i] = 0;
  for (unsigned i = k; i < 3; ++i)
    v[i] = a[i][k];
  v[k] += sign(a[k][k]) * cnorm;
  double rnorm = 0;
  for (unsigned i = k; i < 3; ++i)
    rnorm += square(v[i]);
  rnorm = sqrt(rnorm);
  for (unsigned i = k; i < 3; ++i)
    v[i] /= rnorm;
  return 1;
}

static void reflect_columns(double v[3], double a[3][3], unsigned k)
{
  for (unsigned j = k; j < 3; ++j) {
    double dot = 0;
    for (unsigned i = k; i < 3; ++i)
      dot += a[i][j] * v[i];
    for (unsigned i = k; i < 3; ++i)
      a[i][j] -= 2 * dot * v[i];
  }
}

static void reflect_rows(double v[3], double q[3][3], unsigned k)
{
  for (unsigned i = 0; i < 3; ++i) {
    double dot = 0;
    for (unsigned j = k; j < 3; ++j)
      dot += q[i][j] * v[j];
    for (unsigned j = k; j < 3; ++j)
      q[i][j] -= 2 * dot * v[j];
  }
}

static void qr_step(double a[3][3], unsigned k, double q[3][3])
{
  double v[3];
  if (!get_reflector(a, v, k))
    return;
  reflect_columns(v, a, k);
  reflect_rows(v, q, k);
}

static void copy(double a[3][3], double b[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    b[i][j] = a[i][j];
}

static void fill_identity(double q[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    q[i][j] = ((double)(i==j));
}

void qr_decomp(double a[3][3], double q[3][3], double r[3][3])
{
  copy(a, r);
  fill_identity(q);
  for (unsigned k = 0; k < 3; ++k)
    qr_step(r, k, q);
}
