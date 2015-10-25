#include "shuffle.h"

#include <stdio.h>

#include "comm.h"
#include "global.h"
#include "ints.h"
#include "loop.h"

struct shuffle {
  unsigned nsends;
  unsigned nrecvs;
  unsigned nsend_ents;
  unsigned nrecv_ents;
  unsigned* send_peers;
  unsigned* recv_peers;
  unsigned* send_counts;
  unsigned* recv_counts;
  unsigned* offset_of_sends;
  unsigned* offset_of_recvs;
  unsigned* send_of_ents;
  unsigned* send_idx_of_ents;
  struct comm* c;
};

struct shuffle* new_shuffle(unsigned n, unsigned const* parts)
{
  struct shuffle* s = LOOP_HOST_MALLOC(struct shuffle, 1);
  categorize_by_part(parts, n,
      &s->send_of_ents, &s->send_idx_of_ents,
      &s->nsends, &s->send_peers, &s->send_counts);
  s->c = comm_graph(comm_using(), s->nsends, s->send_peers, s->send_counts);
  loop_free(s->send_peers);
  loop_free(s->send_counts);
  comm_adjacent(s->c,
      &s->nrecvs, &s->recv_peers, &s->recv_counts,
      &s->nsends, &s->send_peers, &s->send_counts);
  s->offset_of_sends = uints_exscan(s->send_counts, s->nsends);
  s->offset_of_recvs = uints_exscan(s->recv_counts, s->nrecvs);
  s->nsend_ents = s->offset_of_sends[s->nsends];
  s->nrecv_ents = s->offset_of_recvs[s->nrecvs];
  return s;
}

void print_shuffle(struct shuffle* s)
{
  printf("shuffle: %u sends, %u receives\n",
      s->nsends, s->nrecvs);
  printf("sends:\n");
  for (unsigned i = 0; i < s->nsends; ++i)
    printf("%u starting at %u to %u\n",
        s->offset_of_sends[i + 1] - s->offset_of_sends[i],
        s->offset_of_sends[i], s->send_peers[i]);
  printf("receives:\n");
  for (unsigned i = 0; i < s->nrecvs; ++i)
    printf("%u starting at %u from %u\n",
        s->offset_of_recvs[i + 1] - s->offset_of_recvs[i],
        s->offset_of_recvs[i], s->recv_peers[i]);
}

void free_shuffle(struct shuffle* s)
{
  loop_free(s->send_peers);
  loop_free(s->recv_peers);
  loop_free(s->send_counts);
  loop_free(s->recv_counts);
  loop_free(s->offset_of_sends);
  loop_free(s->offset_of_recvs);
  loop_free(s->send_of_ents);
  loop_free(s->send_idx_of_ents);
  comm_free(s->c);
  loop_host_free(s);
}

unsigned* shuffle_uints(struct shuffle* s, unsigned const* a, unsigned width)
{
  unsigned* sendbuf = sort_uints_by_category(
      a, width, s->nsend_ents, s->send_of_ents, s->send_idx_of_ents,
      s->offset_of_sends);
  unsigned* recvbuf = LOOP_MALLOC(unsigned, s->nrecv_ents);
  comm_exch_uints(s->c, width, sendbuf, s->send_counts, s->offset_of_sends,
      recvbuf, s->recv_counts, s->offset_of_recvs);
  loop_free(sendbuf);
  return recvbuf;
}

unsigned shuffle_recv_size(struct shuffle* s)
{
  return s->nrecv_ents;
}
