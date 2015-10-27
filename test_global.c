#include <assert.h>
#include <stdio.h>

#include "comm.h"
#include "loop.h"
#include "owners_from_global.h"

int main()
{
  comm_init();
  unsigned* own_parts;
  unsigned* own_idxs;
  unsigned n;
  assert(comm_size() == 2);
  if (comm_rank() == 0) {
    n = 1;
    unsigned long global[1] = {0};
    owners_from_global(n, global, &own_parts, &own_idxs);
  } else {
    n = 1;
    unsigned long global[1] = {1};
    owners_from_global(n, global, &own_parts, &own_idxs);
  }
  for (unsigned i = 0; i < n; ++i)
    printf("%u %u\n", own_parts[i], own_idxs[i]);
  loop_free(own_parts);
  loop_free(own_idxs);
  comm_fini();
}
