#include "inertia.h"

#include <math.h>

#include "algebra.h"
#include "qr.h"

static void cross_matrix(double b[3], double B[3][3])
{
  B[0][0] =     0; B[0][1] = -b[2]; B[0][2] =  b[1];
  B[1][0] =  b[2]; B[1][1] =     0; B[1][2] = -b[0];
  B[2][0] = -b[1]; B[2][1] =  b[0]; B[2][2] =    0;
}

void inertial_contribution(double m, double x[3], double c[3], double ic[3][3])
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
