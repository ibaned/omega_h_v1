#include "inertia.h"
#include "qr.h"
#include "algebra.h"
#include <math.h>
#include <stdio.h>

int main()
{
  unsigned L = 99, W = 100, H = 101;
  printf("L %u W %u H %u\n",L,W,H);
  double c[3] = {
   ((double)L)/2.0,
   ((double)W)/2.0,
   ((double)H)/2.0};
  printf("c %f %f %f\n",c[0],c[1],c[2]);
  double A[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  for (unsigned i = 0; i <= L; ++i)
  for (unsigned j = 0; j <= W; ++j)
  for (unsigned k = 0; k <= H; ++k) {
    double x[3] = { i, j, k };
    double ic[3][3];
    inertial_contribution(1.0, x, c, ic);
    for (unsigned a = 0; a < 3; ++a)
    for (unsigned b = 0; b < 3; ++b)
      A[a][b] += ic[a][b];
  }
  printf("A %f %f %f\n"
         "  %f %f %f\n"
         "  %f %f %f\n",
         A[0][0], A[0][1], A[0][2],
         A[1][0], A[1][1], A[1][2],
         A[2][0], A[2][1], A[2][2]);
  double q[3][3];
  double l[3][3];
  qr_eigen(A, q, l);
  printf("eigenvals:");
  for (unsigned i = 0; i < 3; ++i)
    printf(" %f", l[i][i]);
  printf("\n");
  printf("eigenvectors:\n");
  for (unsigned i = 0; i < 3; ++i) {
    printf("%f %f %f\n",
        q[0][i],
        q[1][i],
        q[2][i]);
  }
  double axis[3];
  least_inertial_axis(A, axis);
  printf("least inertial axis %f %f %f\n",
      axis[0], axis[1], axis[2]);
}
