#include "shuffle.h"

#include <stdio.h>

#include "comm.h"
#include "global.h"
#include "loop.h"

struct shuffle {
  unsigned nsends;
  unsigned nrecvs;
  unsigned* send_peers;
  unsigned* recv_peers;
  unsigned* send_counts;
  unsigned* recv_counts;
  unsigned* ent_peers;
  unsigned* ent_indices;
  struct comm* c;
};

struct shuffle* new_shuffle(unsigned n,
    unsigned const* parts, unsigned const* indices)
{
  (void) indices;
  struct shuffle* s = LOOP_HOST_MALLOC(struct shuffle, 1);
  categorize_by_part(parts, n,
      &s->ent_peers, &s->ent_indices,
      &s->nsends, &s->send_peers, &s->send_counts);
  s->c = comm_graph(comm_using(), s->nsends, s->send_peers, s->send_counts);
  loop_free(s->send_peers);
  loop_free(s->send_counts);
  comm_adjacent(s->c,
      &s->nrecvs, &s->recv_peers, &s->recv_counts,
      &s->nsends, &s->send_peers, &s->send_counts);
  return s;
}

void print_shuffle(struct shuffle* s)
{
  printf("shuffle: %u sends, %u receives\n",
      s->nsends, s->nrecvs);
  printf("sends:\n");
  for (unsigned i = 0; i < s->nsends; ++i)
    printf("%u to %u\n", s->send_counts[i], s->send_peers[i]);
  printf("receives:\n");
  for (unsigned i = 0; i < s->nrecvs; ++i)
    printf("%u from %u\n", s->recv_counts[i], s->recv_peers[i]);
}

void free_shuffle(struct shuffle* s)
{
  loop_free(s->send_peers);
  loop_free(s->recv_peers);
  loop_free(s->send_counts);
  loop_free(s->recv_counts);
  loop_free(s->ent_peers);
  loop_free(s->ent_indices);
  comm_free(s->c);
  loop_host_free(s);
}
