#include <stdio.h>

#include "ints.h"
#include "loop.h"

int main()
{
  unsigned const a[10] = {5,5,4,3,2,1,1,2,3,4};
  unsigned nunique;
  unsigned* unique;
  uints_unique(a, 10, &nunique, &unique);
  for (unsigned i = 0; i < nunique; ++i)
    printf("%u\n", unique[i]);
  loop_free(unique);
}
