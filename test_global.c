#include <stdio.h>

#include "comm.h"
#include "owners_from_global.h"

int main()
{
  comm_init();
  unsigned n = 3;
  unsigned long global[3] = {0,1,2};
  unsigned* own_parts;
  unsigned* own_idxs;
  owners_from_global(n, global, &own_parts, &own_idxs);
  for (unsigned i = 0; i < n; ++i)
    printf("%u %u\n", own_parts[i], own_idxs[i]);
  comm_fini();
}
