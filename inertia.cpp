#include "inertia.hpp"

#include <math.h>

#include "algebra.hpp"
#include "comm.hpp"
#include "doubles.hpp"
#include "loop.hpp"
#include "qr.hpp"

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

/* this function ensures that for the same axis
   we use the same vector (two are possible,
   negatives of one another).
   this helps RIB results look intuitively nice,
   i.e. if we are in 1D then
      1 0 3 2
   is technically a correct answer, but we'd like
      0 1 2 3 */

static void positivize_axis(double a[3])
{
  unsigned signbits = 0;
  for (unsigned i = 0; i < 3; ++i)
    if (a[i] > 0)
      signbits |= (((unsigned)1) << (3-i-1));
  unsigned opp_signbits = signbits ^ 0x7;
  if (opp_signbits > signbits)
    for (unsigned i = 0; i < 3; ++i)
      a[i] = -a[i];
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
  positivize_axis(a);
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

static void find_median_radius(
    unsigned n,
    double const* radii,
    double const* masses,
    double total_mass,
    unsigned is_global,
    unsigned** p_in,
    double* p_wi)
{
  double r = 0;
  double dr = doubles_max(radii, n) / 2;
  if (is_global)
    dr = comm_max_double(dr);
  double hm = total_mass / 2;
  unsigned const fraction_bits = 52;
  for (unsigned i = 0; 1; ++i) {
    unsigned* in = mark_in(n, radii, r);
    double wi = get_weighted_in_count(n, in, masses, is_global);
    if (i == fraction_bits || wi == hm) {
      *p_wi = wi;
      *p_in = in;
      return;
    } else
      loop_free(in);
    if (wi > hm)
      r += dr;
    else
      r -= dr;
    dr /= 2;
  }
}

/* some types of input data have points that happen
   to all be at the same radius, and in the worst
   case they are at the center, meaning that
   unless we change the axis we can't find a good
   median radius.
   this algorithm tries some silly perturbations
   if the original axis can't satisfy the
   imbalance ceiling given */

static void find_median_radius_perturbed(
    unsigned n,
    double const* coords,
    double const* c,
    double const* a,
    double const* masses,
    double total_mass,
    unsigned is_global,
    unsigned** p_in)
{
  /* magic constants... */
  double const max_imb = 0.05;
  double const epsilon = 1e-6;
  for (unsigned i = 0; 1; ++i) {
    double pa[3];
    for (unsigned j = 0; j < 3; ++j)
      pa[j] = a[j] + ((i>>j)&1) * epsilon;
    /* pa[j]=a[j] when i=0 */
    double wi;
    double* radii = get_radii(n, coords, c, pa);
    find_median_radius(n, radii, masses, total_mass, is_global, p_in, &wi);
    loop_free(radii);
    double imb = fabs(2*(wi / total_mass) - 1);
    if (imb <= max_imb || i == 7)
      return;
    else
      loop_free(*p_in);
  }
}

unsigned* mark_inertial_bisection(
    unsigned n,
    double const* coords,
    double const* masses,
    unsigned is_global)
{
  double total_mass;
  if (masses)
    total_mass = doubles_sum(masses, n);
  else
    total_mass = n;
  if (is_global)
    total_mass = comm_add_double(total_mass);
  double c[3];
  get_center_of_mass(n, coords, masses, total_mass, c, is_global);
  double a[3];
  get_axis(n, coords, masses, c, a, is_global);
  unsigned* in;
  find_median_radius_perturbed(n, coords, c, a,
      masses, total_mass, is_global, &in);
  return in;
}
