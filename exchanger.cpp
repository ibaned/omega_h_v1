#include "exchanger.hpp"

#include <cassert>
#include <cstring>
#include <utility> //for std::swap

#include "arrays.hpp"
#include "comm.hpp"
#include "global.hpp"
#include "ints.hpp"
#include "invert_map.hpp"
#include "loop.hpp"
#include "tag.hpp"

namespace omega_h {

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
    unsigned** p_send_order,
    unsigned* p_nsends,
    unsigned** p_send_ranks,
    unsigned** p_send_offsets)
{
  /* queued[i]==1 iff (i) is not part of a message yet */
  unsigned* queued = LOOP_MALLOC(unsigned, nsent);
  for (unsigned i = 0; i < nsent; ++i)
    queued[i] = 1;
  unsigned* send_of_sent = LOOP_MALLOC(unsigned, nsent);
  unsigned* send_order = LOOP_MALLOC(unsigned, nsent);
  unsigned* send_offsets = LOOP_MALLOC(unsigned, nsent + 1);
  unsigned* send_ranks = LOOP_MALLOC(unsigned, nsent);
  /* loop over messages, we don't know how many but
     there certainly won't be more than (nsent) */
  send_offsets[0] = 0;
  unsigned send;
  for (send = 0; send < nsent; ++send) {
    unsigned current_rank = 0;
    unsigned* queue_offsets = uints_exscan(queued, nsent);
    unsigned nqueued = array_at(queue_offsets, nsent);
    if (nqueued == 0) {
      loop_free(queue_offsets);
      break; /* stop when all entries are part of a message */
    }
    /* process the rank of the first queued entry */
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
    send_offsets[send + 1] = send_offsets[send] + array_at(send_idxs, nsent);
    for (unsigned i = 0; i < nsent; ++i)
      if (to_rank[i])
        send_order[i] = send_idxs[i] + send_offsets[send];
    loop_free(to_rank);
    loop_free(send_idxs);
  }
  unsigned nsends = send;
  loop_free(queued);
  *p_send_of_sent = send_of_sent;
  *p_send_order = send_order;
  *p_nsends = nsends;
  /* shrink these arrays to fit, they were
     allocated to the maximum possible size of (nsent) */
  *p_send_ranks = copy_array(send_ranks, nsends);
  loop_free(send_ranks);
  *p_send_offsets = copy_array(send_offsets, nsends + 1);
  loop_free(send_offsets);
}

/* given the number of items to receive,
   and how they are organized to received messages,
   make an array which for every item received lists
   the message it came from.
   this is a slightly embarrassing thing to do, but
   is useful in writing user algorithms */

static unsigned* make_recv_of_recvd(
    unsigned nrecvd,
    unsigned nrecvs,
    unsigned const* recv_offsets)
{
  unsigned* recv_of_recvd = LOOP_MALLOC(unsigned, nrecvd);
  for (unsigned i = 0; i < nrecvs; ++i) {
    unsigned first = recv_offsets[i];
    unsigned end = recv_offsets[i + 1];
    for (unsigned j = first; j < end; ++j)
      recv_of_recvd[j] = i;
  }
  return recv_of_recvd;
}

#define F EX_FOR
#define R EX_REV

struct exchanger* new_exchanger(
    unsigned nsent,
    unsigned const* dest_rank_of_sent)
{
  struct exchanger* ex = LOOP_HOST_MALLOC(struct exchanger, 1);
  memset(ex, 0, sizeof(struct exchanger));
  ex->nitems[F] = nsent;
  sends_from_dest_ranks(nsent, dest_rank_of_sent,
      &ex->msg_of_items[F], &ex->orders[F], &ex->nmsgs[F], &ex->ranks[F],
      &ex->msg_offsets[F]);
  ex->msg_counts[F] = uints_unscan(ex->msg_offsets[F], ex->nmsgs[F]);
  ex->comms[F] = comm_graph(comm_using(), ex->nmsgs[F], ex->ranks[F],
      ex->msg_counts[F]);
  comm_recvs(ex->comms[F], &ex->nmsgs[R], &ex->ranks[R], &ex->msg_counts[R]);
  ex->msg_offsets[R] = uints_exscan(ex->msg_counts[R], ex->nmsgs[R]);
  ex->nitems[R] = ex->msg_offsets[R][ex->nmsgs[R]];
  ex->comms[R] = comm_graph_exact(comm_using(),
      ex->nmsgs[F], ex->ranks[F], ex->msg_counts[F],
      ex->nmsgs[R], ex->ranks[R], ex->msg_counts[R]);
  ex->msg_of_items[R] = make_recv_of_recvd(ex->nitems[R], ex->nmsgs[R],
      ex->msg_offsets[R]);
  return ex;
}

void set_exchanger_dests(
    struct exchanger* ex,
    unsigned ndests,
    unsigned const* dest_idx_of_sent)
{
  unsigned* dest_idx_of_recvd = exchange(ex, 1, dest_idx_of_sent,
      EX_FOR, EX_ITEM);
  ex->nroots[R] = ndests;
  invert_map(ex->nitems[R], dest_idx_of_recvd, ndests,
      &ex->orders[R], &ex->items_of_roots_offsets[R]);
  loop_free(dest_idx_of_recvd);
  unsigned* msg_of_items = reorder_array_inv(ex->msg_of_items[R],
      ex->orders[R], ex->nitems[R], 1);
  loop_free(ex->msg_of_items[R]);
  ex->msg_of_items[R] = msg_of_items;
}

void set_exchanger_srcs(
    struct exchanger* ex,
    unsigned nsrcs,
    unsigned const* sent_of_srcs_offsets)
{
  ex->nroots[F] = nsrcs;
  ex->items_of_roots_offsets[F] = copy_array(sent_of_srcs_offsets, nsrcs + 1);
}

static enum exch_dir opp_dir(enum exch_dir d)
{
  if (d == EX_FOR)
    return EX_REV;
  return EX_FOR;
}

template <typename T>
T* exchange(struct exchanger* ex, unsigned width,
    T const* data, enum exch_dir dir, enum exch_start start)
{
  enum exch_dir odir = opp_dir(dir);
  T const* current = data;
  T* last = 0;
  if (start == EX_ROOT && ex->items_of_roots_offsets[dir]) {
    T* expanded = expand_array(current,
        ex->items_of_roots_offsets[dir], ex->nroots[dir], width);
    current = last = expanded;
  }
  if (ex->orders[dir]) {
    T* orderd = reorder_array(current, ex->orders[dir],
        ex->nitems[dir], width);
    loop_free(last);
    current = last = orderd;
  }
  T* recvd = LOOP_MALLOC(T, ex->nitems[odir] * width);
  comm_exchange<T>(ex->comms[dir], width,
      current, ex->msg_counts[dir], ex->msg_offsets[dir],
      recvd,  ex->msg_counts[odir], ex->msg_offsets[odir]);
  loop_free(last);
  current = last = recvd;
  if (ex->orders[odir]) {
    T* unorderd = reorder_array_inv(current, ex->orders[odir],
        ex->nitems[odir], width);
    loop_free(last);
    current = last = unorderd;
  }
  return last;
}

template unsigned* exchange(struct exchanger* ex, unsigned width,
    unsigned const* data, enum exch_dir dir, enum exch_start start);
template unsigned long* exchange(struct exchanger* ex, unsigned width,
    unsigned long const* data, enum exch_dir dir, enum exch_start start);
template double* exchange(struct exchanger* ex, unsigned width,
    double const* data, enum exch_dir dir, enum exch_start start);

void free_exchanger(struct exchanger* ex)
{
  for (unsigned i = 0; i < 2; ++i) {
    comm_free(ex->comms[i]);
    loop_free(ex->ranks[i]);
    loop_free(ex->msg_counts[i]);
    loop_free(ex->msg_offsets[i]);
    loop_free(ex->orders[i]);
    loop_free(ex->msg_of_items[i]);
    loop_free(ex->items_of_roots_offsets[i]);
  }
  loop_host_free(ex);
}

template <typename T>
static void ex_swap(T a[2])
{
  std::swap((a)[F], (a)[R]);
}

void reverse_exchanger(struct exchanger* ex)
{
  ex_swap(ex->comms);
  ex_swap(ex->nitems);
  ex_swap(ex->nroots);
  ex_swap(ex->nmsgs);
  ex_swap(ex->ranks);
  ex_swap(ex->msg_counts);
  ex_swap(ex->msg_offsets);
  ex_swap(ex->orders);
  ex_swap(ex->msg_of_items);
  ex_swap(ex->items_of_roots_offsets);
}

struct exchanger* make_reverse_exchanger(unsigned nsent, unsigned nrecvd,
    unsigned const* recvd_ranks, unsigned const* recvd_ids)
{
  struct exchanger* ex = new_exchanger(nrecvd, recvd_ranks);
  set_exchanger_dests(ex, nsent, recvd_ids);
  reverse_exchanger(ex);
  return ex;
}

double* exchange_doubles_max(struct exchanger* ex, unsigned width,
    double const* data, enum exch_dir dir, enum exch_start start)
{
  double* to_reduce = exchange(ex, width, data, dir, start);
  enum exch_dir od = opp_dir(dir);
  double* out = LOOP_MALLOC(double, width * ex->nroots[od]);
  doubles_max_into(ex->nroots[od], width, to_reduce,
      ex->items_of_roots_offsets[od], out);
  loop_free(to_reduce);
  return out;
}

}
