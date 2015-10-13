#include "inertia.h"

#include <math.h>

#include "algebra.h"
#include "doubles.h"
#include "loop.h"
#include "qr.h"

static inline void zero_3x3(double a[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    a[i][j] = 0;
}

static inline void add_3x3(double a[3][3], double b[3][3], double c[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    c[i][j] = a[i][j] + b[i][j];
}

static void cross_matrix(double b[3], double B[3][3])
{
  B[0][0] =     0; B[0][1] = -b[2]; B[0][2] =  b[1];
  B[1][0] =  b[2]; B[1][1] =     0; B[1][2] = -b[0];
  B[2][0] = -b[1]; B[2][1] =  b[0]; B[2][2] =    0;
}

void inertia_contribution(
    double m,
    double const* x,
    double const* c,
    double ic[3][3])
{
  double dx[3];
  subtract_vectors(x, c, dx, 3);
  double B[3][3];
  cross_matrix(dx, B);
  mul_3x3(B, B, ic);
  scale_3x3(ic, -m, ic);
}

void least_inertial_axis(double IC[3][3], double a[3])
{
  double q[3][3];
  double l[3][3];
  qr_eigen(IC, q, l);
  unsigned best = 0;
  for (unsigned i = 1; i < 3; ++i)
    if (fabs(l[i][i]) < fabs(l[best][best]))
      best = i;
  for (unsigned i = 0; i < 3; ++i)
    a[i] = q[i][best];
}

static void local_weighted_coords(
    unsigned n,
    double const* coords,
    double const* masses,
    double* c)
{
  for (unsigned i = 0; i < 3; ++i)
    c[i] = 0;
  for (unsigned i = 0; i < n; ++i)
  for (unsigned j = 0; j < 3; ++j) {
    double m = masses ? masses[i] : 1;
    c[j] += coords[i * 3 + j] * m;
  }
}

static void local_center_of_mass(
    unsigned n,
    double const* coords,
    double const* masses,
    double lm,
    double* c)
{
  local_weighted_coords(n, coords, masses, c);
  for (unsigned i = 0; i < 3; ++i)
    c[i] /= lm;
}

static void local_inertia(
    unsigned n,
    double const* coords,
    double const* masses,
    double const* c,
    double ic[3][3])
{
  zero_3x3(ic);
  for (unsigned i = 0; i < n; ++i) {
    double pic[3][3];
    double m = masses ? masses[i] : 1;
    inertia_contribution(m, coords + i * 3, c, pic);
    add_3x3(ic, pic, ic);
  }
}

static void local_axis(
    unsigned n,
    double const* coords,
    double const* masses,
    double const* c,
    double* a)
{
  double ic[3][3];
  local_inertia(n, coords, masses, c, ic);
  least_inertial_axis(ic, a);
}

static double* local_radii(
    unsigned n,
    double const* coords,
    double const* c,
    double const* a)
{
  double* out = LOOP_MALLOC(double, n);
  for (unsigned i = 0; i < n; ++i) {
    double x[3];
    subtract_vectors(coords + i * 3, c, x, 3);
    double pr = dot_product(x, a, 3);
    out[i] = pr;
  }
  return out;
}

static unsigned* mark_local_in(
    unsigned n,
    double const* radii,
    double r)
{
  unsigned* in = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i) {
    in[i] = (radii[i] >= r);
  }
  return in;
}

static double local_weighted_in(
    unsigned n,
    unsigned const* in,
    double const* masses)
{
  double s = 0;
  for (unsigned i = 0; i < n; ++i)
    if (in[i]) {
      double m = masses ? masses[i] : 1;
      s += m;
    }
  return s;
}

static double local_mean_radius(
    unsigned n,
    double const* radii,
    double const* masses,
    double local_mass)
{
  double r = 0;
  double dr = doubles_max(radii, n) / 2;
  double hm = local_mass / 2;
  for (unsigned i = 0; i < 52; ++i) {
    unsigned* in = mark_local_in(n, radii, r);
    double wi = local_weighted_in(n, in, masses);
    LOOP_FREE(in);
    if (wi > hm)
      r += dr;
    else
      r -= dr;
    dr /= 2;
  }
  return r;
}

unsigned* local_inertia_mark(
    unsigned n,
    double const* coords,
    double const* masses)
{
  double lm;
  if (masses)
    lm = doubles_sum(masses, n);
  else
    lm = n;
  double c[3];
  local_center_of_mass(n, coords, masses, lm, c);
  double a[3];
  local_axis(n, coords, masses, c, a);
  double* radii = local_radii(n, coords, c, a);
  double r = local_mean_radius(n, radii, masses, lm);
  unsigned* in = mark_local_in(n, radii, r);
  LOOP_FREE(radii);
  return in;
}
