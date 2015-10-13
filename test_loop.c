#include "loop.h"

LOOP_KERNEL(fill, double* a)
  a[i] = 0;
}

int main()
{
  unsigned n = 1024*1024;
  double* a = LOOP_MALLOC(double, n);
  LOOP_EXEC(fill, n, a);
  LOOP_FREE(a);
}
