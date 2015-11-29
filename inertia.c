#include "inertia.h"

#include <math.h>

#include "algebra.h"
#include "comm.h"
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

static void get_weighted_coords(
    unsigned n,
    double const* coords,
    double const* masses,
    double* c,
    unsigned is_global)
{
  for (unsigned i = 0; i < 3; ++i)
    c[i] = 0;
  for (unsigned i = 0; i < n; ++i)
  for (unsigned j = 0; j < 3; ++j) {
    double m = masses ? masses[i] : 1;
    c[j] += coords[i * 3 + j] * m;
  }
  if (is_global)
    comm_add_doubles(c, 3);
}

static void get_center_of_mass(
    unsigned n,
    double const* coords,
    double const* masses,
    double lm,
    double* c,
    unsigned is_global)
{
  get_weighted_coords(n, coords, masses, c, is_global);
  for (unsigned i = 0; i < 3; ++i)
    c[i] /= lm;
}

static void get_total_inertia(
    unsigned n,
    double const* coords,
    double const* masses,
    double const* c,
    double ic[3][3],
    unsigned is_global)
{
  zero_3x3(ic);
  for (unsigned i = 0; i < n; ++i) {
    double pic[3][3];
    double m = masses ? masses[i] : 1;
    inertia_contribution(m, coords + i * 3, c, pic);
    add_3x3(ic, pic, ic);
  }
  if (is_global)
    comm_add_doubles(ic[0], 9);
}

static void get_axis(
    unsigned n,
    double const* coords,
    double const* masses,
    double const* c,
    double* a,
    unsigned is_global)
{
  double ic[3][3];
  get_total_inertia(n, coords, masses, c, ic, is_global);
  least_inertial_axis(ic, a);
}

static double* get_radii(
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

static unsigned* mark_in(
    unsigned n,
    double const* radii,
    double r)
{
  unsigned* in = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    in[i] = (radii[i] >= r);
  return in;
}

static double get_weighted_in_count(
    unsigned n,
    unsigned const* in,
    double const* masses,
    unsigned is_global)
{
  double s = 0;
  for (unsigned i = 0; i < n; ++i)
    if (in[i]) {
      double m = masses ? masses[i] : 1;
      s += m;
    }
  if (is_global)
    return comm_add_double(s);
  return s;
}

static double get_mean_radius(
    unsigned n,
    double const* radii,
    double const* masses,
    double total_mass,
    unsigned is_global)
{
  double r = 0;
  double dr = doubles_max(radii, n) / 2;
  double hm = total_mass / 2;
  unsigned const fraction_bits = 52;
  for (unsigned i = 0; i < fraction_bits; ++i) {
    unsigned* in = mark_in(n, radii, r);
    double wi = get_weighted_in_count(n, in, masses, is_global);
    loop_free(in);
    if (wi > hm)
      r += dr;
    else
      r -= dr;
    dr /= 2;
  }
  return r;
}

unsigned* mark_inertial_bisection(
    unsigned n,
    double const* coords,
    double const* masses,
    unsigned is_global)
{
  double tm;
  if (masses)
    tm = doubles_sum(masses, n);
  else
    tm = n;
  if (is_global)
    tm = comm_add_double(tm);
  double c[3];
  get_center_of_mass(n, coords, masses, tm, c, is_global);
  double a[3];
  get_axis(n, coords, masses, c, a, is_global);
  double* radii = get_radii(n, coords, c, a);
  double r = get_mean_radius(n, radii, masses, tm, is_global);
  unsigned* in = mark_in(n, radii, r);
  loop_free(radii);
  return in;
}
