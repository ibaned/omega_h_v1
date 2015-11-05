#include "exchanger.h"

#include "comm.h"
#include "global.h"
#include "ints.h"
#include "loop.h"

struct exchanger* new_exchanger(unsigned nsent,
    unsigned const* dest_rank_of_sent)
{
  struct exchanger* ex = LOOP_HOST_MALLOC(struct exchanger, 1);
  ex->nsent = nsent;
  categorize_by_part(dest_rank_of_sent, nsent,
      &ex->send_of_sent, &ex->send_idx_of_sent,
      &ex->nsends, &ex->send_ranks, &ex->send_counts);
  ex->forward_comm = comm_graph(comm_using(), ex->nsends, ex->send_ranks,
      ex->send_counts);
  comm_recvs(ex->forward_comm, &ex->nrecvs, &ex->recv_ranks, &ex->recv_counts);
  ex->send_offsets = uints_exscan(ex->send_counts, ex->nsends);
  ex->recv_offsets = uints_exscan(ex->recv_counts, ex->nrecvs);
  ex->nrecvd = ex->recv_offsets[ex->nrecvs];
  ex->reverse_comm = comm_graph_exact(comm_using(),
      ex->nsends, ex->send_ranks, ex->send_counts,
      ex->nrecvs, ex->recv_ranks, ex->recv_counts);
  return ex;
}

unsigned* exchange_uints(struct exchanger* ex, unsigned width,
    unsigned const* sent)
{
  unsigned* sorted = sort_uints_by_category(sent, width, ex->nsent,
      ex->send_of_sent, ex->send_idx_of_sent, ex->send_offsets);
  unsigned* recvd = LOOP_HOST_MALLOC(unsigned, ex->nrecvd * width);
  comm_exch_uints(ex->forward_comm, width,
      sorted, ex->send_counts, ex->send_offsets,
      recvd,  ex->recv_counts, ex->recv_offsets);
  loop_free(sorted);
  return recvd;
}

unsigned* unexchange_uints(struct exchanger* ex, unsigned width,
    unsigned const* recvd)
{
  unsigned* sorted = LOOP_HOST_MALLOC(unsigned, ex->nsent * width);
  comm_exch_uints(ex->reverse_comm, width,
      recvd,  ex->recv_counts, ex->recv_offsets,
      sorted, ex->send_counts, ex->send_offsets);
  unsigned* sent = unsort_uints_by_category(sorted, width, ex->nsent,
      ex->send_of_sent, ex->send_idx_of_sent, ex->send_offsets);
  loop_free(sorted);
  return sent;
}

void free_exchanger(struct exchanger* ex)
{
  comm_free(ex->forward_comm);
  comm_free(ex->reverse_comm);
  loop_free(ex->send_ranks);
  loop_free(ex->recv_ranks);
  loop_free(ex->send_counts);
  loop_free(ex->recv_counts);
  loop_free(ex->send_offsets);
  loop_free(ex->recv_offsets);
  loop_free(ex->send_of_sent);
  loop_free(ex->send_idx_of_sent);
  loop_host_free(ex);
}
