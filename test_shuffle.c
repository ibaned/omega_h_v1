#include "shuffle.h"

#include <assert.h>

#include "comm.h"

int main()
{
  comm_init();
  assert(comm_size() == 2);
  struct shuffle* s;
  if (comm_rank() == 0) {
    unsigned n = 1;
    unsigned const parts[1] = {1};
    unsigned const indices[1] = {0};
    s = new_shuffle(n, parts, indices);
  } else {
    unsigned n = 1;
    unsigned const parts[1] = {0};
    unsigned const indices[1] = {0};
    s = new_shuffle(n, parts, indices);
  }
  print_shuffle(s);
  free_shuffle(s);
  comm_fini();
}
