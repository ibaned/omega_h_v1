#include "inertia.h"
#include "algebra.h"
#include "find_roots.h"

static void cross_matrix(double b[3], double B[3][3])
{
  B[0][0] =     0; B[0][1] = -b[2]; B[0][2] =  b[1];
  B[1][0] =  b[2]; B[1][1] =     0; B[1][2] = -b[0];
  B[2][0] = -b[1]; B[2][1] =  b[0]; B[2][2] =    0;
}

static void mmm_3x3(double a[3][3], double b[3][3], double c[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j) {
    c[i][j] = 0;
    for (unsigned k = 0; k < 3; ++k)
      c[i][j] += a[i][k] * b[k][j];
  }
}

void inertial_contribution(double m, double x[3], double c[3], double ic[3][3])
{
  double dx[3];
  subtract_vectors(x, c, dx, 3);
  double B[3][3];
  cross_matrix(dx, B);
  mmm_3x3(B, B, ic);
  scale_3x3(ic, -m);
}

static double trace_3x3(double A[3][3])
{
  return A[0][0] + A[1][1] + A[2][2];
}

static void char_poly_3x3(double A[3][3],
    double* a, double* b, double* c, double* d)
{
  *a = -1;
  double tA = trace_3x3(A);
  *b = tA;
  double Asq[3][3];
  mmm_3x3(A, A, Asq);
  *c = -(tA * tA - trace_3x3(Asq)) / 2.;
  *d = det_3x3(A);
}

unsigned eigenvals_3x3(double A[3][3], double l[])
{
  double a,b,c,d;
  char_poly_3x3(A, &a, &b, &c, &d);
  return find_cubic_roots(a, b, c, d, l);
}
