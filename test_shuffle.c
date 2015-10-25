#include "shuffle.h"

#include <assert.h>
#include <stdio.h>

#include "comm.h"
#include "loop.h"

int main()
{
  comm_init();
  assert(comm_size() == 2);
  if (comm_rank() == 0) {
    unsigned n = 3;
    unsigned const parts[3] = {1,0,1};
    struct shuffle* s = new_shuffle(n, parts);
    print_shuffle(s);
    unsigned const sent[3] = {1,2,3};
    unsigned* recvd = shuffle_uints(s, sent);
    printf("recvd:\n");
    unsigned nr = shuffle_recv_size(s);
    for (unsigned i = 0; i < nr; ++i)
      printf("%u\n", recvd[i]);
    loop_free(recvd);
    free_shuffle(s);
  } else {
    unsigned n = 2;
    unsigned const parts[2] = {0,0};
    struct shuffle* s = new_shuffle(n, parts);
    print_shuffle(s);
    unsigned const sent[2] = {4,5};
    unsigned* recvd = shuffle_uints(s, sent);
    printf("recvd:\n");
    unsigned nr = shuffle_recv_size(s);
    for (unsigned i = 0; i < nr; ++i)
      printf("%u\n", recvd[i]);
    loop_free(recvd);
    free_shuffle(s);
  }
  comm_fini();
}
