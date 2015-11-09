#include <assert.h>
#include <stdio.h>

#include "comm.h"
#include "loop.h"
#include "owners_from_global.h"

static void test_owners_from_global(void)
{
  unsigned* own_parts;
  unsigned* own_idxs;
  unsigned n;
  assert(comm_size() == 2);
  if (comm_rank() == 0) {
    n = 2;
    unsigned long global[2] = {0,1};
    owners_from_global(n, global, &own_parts, &own_idxs);
    assert(own_parts[0] == 0);
    assert(own_parts[1] == 0);
    assert(own_idxs[0] == 0);
    assert(own_idxs[1] == 1);
  } else {
    n = 2;
    unsigned long global[2] = {1,2};
    owners_from_global(n, global, &own_parts, &own_idxs);
    assert(own_parts[0] == 0);
    assert(own_parts[1] == 1);
    assert(own_idxs[0] == 1);
    assert(own_idxs[1] == 1);
  }
  loop_free(own_parts);
  loop_free(own_idxs);
}

int main()
{
  comm_init();
  test_owners_from_global();
  comm_fini();
}
