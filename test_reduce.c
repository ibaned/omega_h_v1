#include "arrays.h"

#include <stdio.h>

int main()
{
  unsigned offsets[3] = {0,2,5};
  double data[5] = {11,22,3,7,1};
  double output[2];
  doubles_max_into(2, 1, data, offsets, output);
  printf("output: {%f %f}\n", output[0], output[1]);
}
