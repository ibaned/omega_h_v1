#include "shuffle.h"

#include <assert.h>

#include "comm.h"

int main()
{
  comm_init();
  assert(comm_size() == 2);
  struct shuffle* s;
  if (comm_rank() == 0) {
    unsigned n = 3;
    unsigned const parts[3] = {1,1,0};
    unsigned const indices[3] = {0,1,2};
    s = new_shuffle(n, parts, indices);
  } else {
    unsigned n = 2;
    unsigned const parts[2] = {0,0};
    unsigned const indices[2] = {0,1};
    s = new_shuffle(n, parts, indices);
  }
  print_shuffle(s);
  free_shuffle(s);
  comm_fini();
}
