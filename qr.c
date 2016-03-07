#include "qr.h"

#include <assert.h>
#include <math.h>

#include "algebra.h"

static double sign(double x)
{
  return (x < 0) ? -1 : 1;
}

static double square(double x)
{
  return x * x;
}

static unsigned get_reflector(double a[3][3], double v[3], unsigned k,
    unsigned o)
{
  double cnorm = 0;
  for (unsigned i = k + o; i < 3; ++i)
    cnorm += square(a[i][k]);
  cnorm = sqrt(cnorm);
  if (cnorm < 1e-14)
    return 0;
  for (unsigned i = 0; i < k + o; ++i)
    v[i] = 0;
  for (unsigned i = k + o; i < 3; ++i)
    v[i] = a[i][k];
  v[k + o] += sign(a[k + o][k]) * cnorm;
  double rnorm = 0;
  for (unsigned i = k + o; i < 3; ++i)
    rnorm += square(v[i]);
  rnorm = sqrt(rnorm);
  for (unsigned i = k + o; i < 3; ++i)
    v[i] /= rnorm;
  return 1;
}

static unsigned get_reflector2(double a[MAX_PTS][4],
    double v[MAX_PTS], unsigned k, unsigned npts)
{
  double cnorm = 0;
  for (unsigned i = k; i < npts; ++i)
    cnorm += square(a[i][k]);
  cnorm = sqrt(cnorm);
  if (cnorm < 1e-14)
    return 0;
  for (unsigned i = 0; i < k; ++i)
    v[i] = 0;
  for (unsigned i = k; i < npts; ++i)
    v[i] = a[i][k];
  v[k] += sign(a[k][k]) * cnorm;
  double rnorm = 0;
  for (unsigned i = k; i < npts; ++i)
    rnorm += square(v[i]);
  rnorm = sqrt(rnorm);
  for (unsigned i = k; i < npts; ++i)
    v[i] /= rnorm;
  return 1;
}

static void reflect_columns(double v[3], double a[3][3], unsigned k,
    unsigned o)
{
  for (unsigned j = 0; j < 3; ++j) {
    double dot = 0;
    for (unsigned i = k + o; i < 3; ++i)
      dot += a[i][j] * v[i];
    for (unsigned i = k + o; i < 3; ++i)
      a[i][j] -= 2 * dot * v[i];
  }
}

static void reflect_columns2(double v[MAX_PTS], double a[MAX_PTS][4],
    unsigned k, unsigned npts)
{
  for (unsigned j = 0; j < 4; ++j) {
    double dot = 0;
    for (unsigned i = k; i < npts; ++i)
      dot += a[i][j] * v[i];
    for (unsigned i = k; i < npts; ++i)
      a[i][j] -= 2 * dot * v[i];
  }
}

static void reflect_rows(double v[3], double q[3][3], unsigned k,
    unsigned o)
{
  for (unsigned i = 0; i < 3; ++i) {
    double dot = 0;
    for (unsigned j = k + o; j < 3; ++j)
      dot += q[i][j] * v[j];
    for (unsigned j = k + o; j < 3; ++j)
      q[i][j] -= 2 * dot * v[j];
  }
}

static void reflect_rows2(double v[MAX_PTS], double q[MAX_PTS][MAX_PTS],
    unsigned k, unsigned npts)
{
  for (unsigned i = 0; i < npts; ++i) {
    double dot = 0;
    for (unsigned j = k; j < npts; ++j)
      dot += q[i][j] * v[j];
    for (unsigned j = k; j < npts; ++j)
      q[i][j] -= 2 * dot * v[j];
  }
}

static void copy(double a[3][3], double b[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    b[i][j] = a[i][j];
}

static void copy2(double a[MAX_PTS][4], double b[MAX_PTS][4], unsigned npts)
{
  for (unsigned i = 0; i < npts; ++i)
  for (unsigned j = 0; j < 4; ++j)
    b[i][j] = a[i][j];
}

static void fill_identity(double q[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    q[i][j] = ((double)(i==j));
}

static void fill_identity2(double q[MAX_PTS][MAX_PTS], unsigned npts)
{
  for (unsigned i = 0; i < npts; ++i)
  for (unsigned j = 0; j < npts; ++j)
    q[i][j] = ((double)(i==j));
}

void qr_decomp(double a[3][3], double q[3][3], double r[3][3])
{
  copy(a, r);
  fill_identity(q);
  double v[3] = {0};
  for (unsigned k = 0; k < 3; ++k)
    if (get_reflector(r, v, k, 0)) {
      reflect_columns(v, r, k, 0);
      reflect_rows(v, q, k, 0);
    }
}

void qr_decomp2(
    double a[MAX_PTS][4],
    double q[MAX_PTS][MAX_PTS],
    double r[MAX_PTS][4],
    unsigned npts)
{
  copy2(a, r, npts);
  fill_identity2(q, npts);
  double v[MAX_PTS] = {0};
  for (unsigned k = 0; k < 4; ++k)
    if (get_reflector2(r, v, k, npts)) {
      reflect_columns2(v, r, k, npts);
      reflect_rows2(v, q, k, npts);
    }
}

static void hessenberg(double a[3][3], double q[3][3], double h[3][3])
{
  copy(a, h);
  fill_identity(q);
  double v[3] = {0};
  if (get_reflector(h, v, 0, 1)) {
    reflect_columns(v, h, 0, 1);
    reflect_rows(v, h, 0, 1);
    reflect_rows(v, q, 0, 1);
  }
}

static unsigned reduce(double a[3][3], unsigned* m)
{
  for (; *m > 1; --(*m))
    if (fabs(a[*m - 2][*m - 1]) > 1e-10 ||
        fabs(a[*m - 1][*m - 2]) > 1e-10)
      return 1;
  return 0;
}

static double wilkinson(double a[3][3], unsigned m)
{
  double amm1 = a[m - 2][m - 2];
  double am = a[m - 1][m - 1];
  double bmm1 = a[m - 2][m - 1];
  double sig = (amm1 - am) / 2;
  double denom = fabs(sig) + sqrt(square(sig) + square(bmm1));
  assert(denom > 1e-10);
  return am - ((sign(sig) * square(bmm1)) / denom);
}

static void shift(double a[3][3], double mu)
{
  for (unsigned i = 0; i < 3; ++i)
    a[i][i] -= mu;
}

void qr_eigen(double a[3][3], double q[3][3], double l[3][3])
{
  hessenberg(a, q, l);
  double rk[3][3];
  double qk[3][3];
  double q2[3][3];
  unsigned m = 3;
  for (unsigned i = 0; i < 100; ++i) {
    if (!reduce(l, &m))
      return;
    double mu = wilkinson(l, m);
    shift(l, mu);
    qr_decomp(l, qk, rk);
    mul_3x3(rk, qk, l);
    shift(l, -mu);
    mul_3x3(q, qk, q2);
    copy(q2, q);
  }
  assert(0);
}
