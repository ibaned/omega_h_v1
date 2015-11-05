#include "exchanger.h"

#include "comm.h"
#include "global.h"
#include "ints.h"
#include "loop.h"

/* given an array that indicates which rank an
   entry is going to,
   this function organizes them into one message
   per destination rank.
   we rely on the assumption that the number of
   messages (nsends) is small to settle on a runtime
   that is O(nsent * nsends * log(nsent)).
   the log(nsent) comes from the runtime of uints_exscan.
*/

static void sends_from_dest_ranks(
    unsigned nsent,
    unsigned const* dest_rank_of_sent,
    unsigned** p_send_of_sent,
    unsigned** p_send_idx_of_sent,
    unsigned* p_nsends,
    unsigned** p_send_ranks,
    unsigned** p_send_counts)
{
  /* queued[i]==1 iff (i) has not been categorized yet */
  unsigned* queued = LOOP_MALLOC(unsigned, nsent);
  for (unsigned i = 0; i < nsent; ++i)
    queued[i] = 1;
  /* loop over categories, we don't know how many but
     there certainly won't be more than (nsent) */
  unsigned* send_of_sent = LOOP_MALLOC(unsigned, nsent);
  unsigned* send_idx_of_sent = LOOP_MALLOC(unsigned, nsent);
  unsigned* send_counts = LOOP_MALLOC(unsigned, nsent);
  unsigned* send_ranks = LOOP_MALLOC(unsigned, nsent);
  unsigned send;
  for (send = 0; send < nsent; ++send) {
    unsigned current_rank = 0;
    unsigned* queue_offsets = uints_exscan(queued, nsent);
    unsigned nqueued = queue_offsets[nsent];
    if (nqueued == 0) {
      loop_free(queue_offsets);
      break; /* stop when all entries categorized */
    }
    /* process the rank of the first uncategorized entry */
    for (unsigned i = 0; i < nsent; ++i)
      if ((queue_offsets[i + 1] - queue_offsets[i] == 1) &&
          queue_offsets[i] == 0)
        current_rank = dest_rank_of_sent[i];
    send_ranks[send] = current_rank;
    loop_free(queue_offsets);
    unsigned* to_rank = LOOP_MALLOC(unsigned, nsent);
    for (unsigned i = 0; i < nsent; ++i) {
      if (dest_rank_of_sent[i] == current_rank) {
        send_of_sent[i] = send;
        to_rank[i] = 1;
        queued[i] = 0;
      } else {
        to_rank[i] = 0;
      }
    }
    unsigned* send_idxs = uints_exscan(to_rank, nsent);
    send_counts[send] = send_idxs[nsent];
    for (unsigned i = 0; i < nsent; ++i)
      if (to_rank[i])
        send_idx_of_sent[i] = send_idxs[i];
    loop_free(to_rank);
    loop_free(send_idxs);
  }
  unsigned nsends = send;
  loop_free(queued);
  *p_send_of_sent = send_of_sent;
  *p_send_idx_of_sent = send_idx_of_sent;
  *p_nsends = nsends;
  /* shrink these arrays to fit, they were
     allocated to the maximum possible size of (nsent) */
  *p_send_ranks = uints_copy(send_ranks, nsends);
  loop_free(send_ranks);
  *p_send_counts = uints_copy(send_counts, nsends);
  loop_free(send_counts);
}

struct exchanger* new_exchanger(unsigned nsent,
    unsigned const* dest_rank_of_sent)
{
  struct exchanger* ex = LOOP_HOST_MALLOC(struct exchanger, 1);
  ex->nsent = nsent;
  sends_from_dest_ranks(nsent, dest_rank_of_sent,
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
