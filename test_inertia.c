#include "inertia.h"
#include "algebra.h"
#include <math.h>
#include <stdio.h>

int main()
{
  unsigned L = 100, W = 100, H = 0;
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
  //printf("x %f %f %f\n",x[0],x[1],x[2]);
    double ic[3][3];
    inertial_contribution(1.0, x, c, ic);
  //printf("ic %f %f %f\n"
  //       "   %f %f %f\n"
  //       "   %f %f %f\n",
  //       ic[0][0], ic[0][1], ic[0][2],
  //       ic[1][0], ic[1][1], ic[1][2],
  //       ic[2][0], ic[2][1], ic[2][2]);
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
  double Amax = inf_norm_3x3(A);
  printf("max A component %f\n", Amax);
  if (Amax != 0.0)
    scale_3x3(A, 1.0 / Amax, A);
  printf("after normalizing:\n"
         "A %f %f %f\n"
         "  %f %f %f\n"
         "  %f %f %f\n",
         A[0][0], A[0][1], A[0][2],
         A[1][0], A[1][1], A[1][2],
         A[2][0], A[2][1], A[2][2]);
  double l[3];
  unsigned nl = eigenvals_3x3(A, l);
  printf("eigenvals:");
  for (unsigned i = 0; i < nl; ++i)
    printf(" %f", l[i]);
  printf("\n");
  printf("eigenvectors:\n");
  for (unsigned i = 0; i < nl; ++i) {
    double v[3];
    eigenvector_3x3(A, l[i], v);
    printf("%f %f %f\n",v[0],v[1],v[2]);
  }
}
