#include <stdio.h>

#include "qr.h"

int main()
{
  double A[3][3] = {
    {0.8147,    0.9134,    0.2785},
    {0.9058,    0.6324,    0.5469},
    {0.1270,    0.0975,    0.9575}
  };
  double q[3][3];
  double l[3][3];
  qr_eigen(A, q, l);
  printf("qr eigenvals:");
  for (unsigned i = 0; i < 3; ++i)
    printf(" %f", l[i][i]);
  printf("\n");
}
