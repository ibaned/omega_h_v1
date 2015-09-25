#include "find_roots.h"
#include <complex.h>
#include <math.h>

unsigned find_linear_root(double a, double b, double* root)
{
  if (a == 0)
    return 0;
  *root = -b / a;
  return 1;
}

unsigned find_quadratic_roots(
    double a, double b, double c,
    double* roots)
{
  if (a == 0)
    return find_linear_root(b, c, roots);
  double D = b * b - 4 * a * c;
  if (D < 0)
    return 0;
  if (D > 0) {
    roots[0] = (-b + sqrt(D)) / (2 * a);
    roots[1] = (-b - sqrt(D)) / (2 * a);
    return 2;
  }
  roots[0] = -b / (2 * a);
  return 1;
}

unsigned find_cubic_roots(
    double a, double b, double c, double d,
    double* roots)
{
  if (a == 0)
    return find_quadratic_roots(b, c, d, roots);
  double D0 = b * b - 3 * a * c;
  double D = 18 * a * b * c * d
           -  4 * b * b * b * d
           +      b * b * c * c
           -  4 * a * c * c * c
           - 27 * a * a * d * d;
  if (D != 0) {
    double D1 =  2 * b * b * b
              -  9 * a * b * c
              + 27 * a * a * d;
    double complex v;
    if (D0 != 0)
      v = csqrt(D1 * D1 - 4 * D0 * D0 * D0);
    else
      v = D1;
    double complex w = (D1 + v) / 2;
    double complex C = cpow(w, 1.0 / 3.0);
    double complex x[3];
    double complex const u[3] = {
      1.0,
      (-1.0 + I * sqrt(3)) / 2.0,
      (-1.0 - I * sqrt(3)) / 2.0
    };
    for (unsigned k = 0; k < 3; ++k)
      x[k] = (-1.0 / (3 * a))
           * (b + u[k] * C + (D0 / (u[k] * C)));
    if (D > 0) {
      for (unsigned k = 0; k < 3; ++k)
        roots[k] = creal(x[k]);
      return 3;
    }
    roots[0] = creal(x[0]);
    for (unsigned k = 1; k < 3; ++k)
      if (fabs(creal(x[k])) > fabs(roots[0]))
        roots[k] = creal(x[k]);
    return 1;
  }
  if (D0 != 0) {
    roots[0] = (9 * a * d - b * c) / (2 * D0);
    roots[1] = ( 4 * a * b * c
               - 9 * a * a * d
               -     b * b * b)
             / (a * D0);
    return 2;
  }
  roots[0] = - b / (3 * a);
  return 1;
}
