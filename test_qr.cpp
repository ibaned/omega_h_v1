#include "qr.hpp"

#include <cstdio>

int main()
{
  double a[MAX_PTS][4] = {
    {1, 0, 0, 0},
    {1, 1, 0, 0},
    {1, 0, 2, 0},
    {1, 0, 0, 3}
  };
  double b[MAX_PTS] = {
    0,
    1,
    2,
    3
  };
  double x[4];
  least_squares_fit(a, b, 4, x);
  for (unsigned i = 0; i < 4; ++i)
    printf("%f\n", x[i]);
}
