#include <stdio.h>

#include "loop.h"
#include "global.h"

int main()
{
  unsigned long global[7] = {0,1,2,3,4,5,6};
  unsigned* part;
  unsigned* local;
  global_to_linpart(global, 7, 7, 4, &part, &local);
  for (unsigned i = 0; i < 7; ++i)
    printf("part %u local %u\n", part[i], local[i]);
  loop_free(part);
  loop_free(local);
}
