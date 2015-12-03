#include "exchanger.h"

#include <assert.h>
#include <string.h>

#include "arrays.h"
#include "comm.h"
#include "global.h"
#include "ints.h"
#include "invert_map.h"
#include "loop.h"
#include "tag.h"

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
    unsigned** p_send_shuffle,
    unsigned* p_nsends,
    unsigned** p_send_ranks,
    unsigned** p_send_offsets)
{
  /* queued[i]==1 iff (i) is not part of a message yet */
  unsigned* queued = LOOP_MALLOC(unsigned, nsent);
  for (unsigned i = 0; i < nsent; ++i)
    queued[i] = 1;
  unsigned* send_of_sent = LOOP_MALLOC(unsigned, nsent);
  unsigned* send_shuffle = LOOP_MALLOC(unsigned, nsent);
  unsigned* send_offsets = LOOP_MALLOC(unsigned, nsent + 1);
  unsigned* send_ranks = LOOP_MALLOC(unsigned, nsent);
  /* loop over messages, we don't know how many but
     there certainly won't be more than (nsent) */
  send_offsets[0] = 0;
  unsigned send;
  for (send = 0; send < nsent; ++send) {
    unsigned current_rank = 0;
    unsigned* queue_offsets = uints_exscan(queued, nsent);
    unsigned nqueued = queue_offsets[nsent];
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
    send_offsets[send + 1] = send_offsets[send] + send_idxs[nsent];
    for (unsigned i = 0; i < nsent; ++i)
      if (to_rank[i])
        send_shuffle[i] = send_idxs[i] + send_offsets[send];
    loop_free(to_rank);
    loop_free(send_idxs);
  }
  unsigned nsends = send;
  loop_free(queued);
  *p_send_of_sent = send_of_sent;
  *p_send_shuffle = send_shuffle;
  *p_nsends = nsends;
  /* shrink these arrays to fit, they were
     allocated to the maximum possible size of (nsent) */
  *p_send_ranks = uints_copy(send_ranks, nsends);
  loop_free(send_ranks);
  *p_send_offsets = uints_copy(send_offsets, nsends + 1);
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
      &ex->msg_of_items[F], &ex->shuffles[F], &ex->nmsgs[F], &ex->ranks[F],
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
  unsigned* dest_idx_of_recvd = exchange_uints(ex, 1, dest_idx_of_sent,
      EX_FOR, EX_ITEM);
  ex->nroots[R] = ndests;
  invert_map(ex->nitems[R], dest_idx_of_recvd, ndests,
      &ex->shuffles[R], &ex->items_of_roots_offsets[R]);
  loop_free(dest_idx_of_recvd);
  unsigned* msg_of_items = uints_unshuffle(ex->nitems[R],
      ex->msg_of_items[R], 1, ex->shuffles[R]);
  loop_free(ex->msg_of_items[R]);
  ex->msg_of_items[R] = msg_of_items;
}

void set_exchanger_srcs(
    struct exchanger* ex,
    unsigned nsrcs,
    unsigned const* sent_of_srcs_offsets)
{
  ex->nroots[F] = nsrcs;
  ex->items_of_roots_offsets[F] = uints_copy(sent_of_srcs_offsets, nsrcs + 1);
}

static enum exch_dir opp_dir(enum exch_dir d)
{
  if (d == EX_FOR)
    return EX_REV;
  return EX_FOR;
}

#define GENERIC_EXCHANGE(T, name) \
T* exchange_##name(struct exchanger* ex, unsigned width, \
    T const* data, enum exch_dir dir, enum exch_start start) \
{ \
  enum exch_dir odir = opp_dir(dir); \
  T const* current = data; \
  T* last = 0; \
  if (start == EX_ROOT) { \
    T* expanded = name##_expand(ex->nroots[dir], current, width, \
        ex->items_of_roots_offsets[dir]); \
    current = last = expanded; \
  } \
  if (ex->shuffles[dir]) { \
    T* shuffled = name##_shuffle(ex->nitems[dir], current, width, \
        ex->shuffles[dir]); \
    loop_free(last); \
    current = last = shuffled; \
  } \
  T* recvd = LOOP_MALLOC(T, ex->nitems[odir] * width); \
  comm_exch_##name(ex->comms[dir], width, \
      current, ex->msg_counts[dir], ex->msg_offsets[dir], \
      recvd,  ex->msg_counts[odir], ex->msg_offsets[odir]); \
  loop_free(last); \
  current = last = recvd; \
  if (ex->shuffles[odir]) { \
    T* unshuffled = name##_unshuffle(ex->nitems[odir], current, width, \
        ex->shuffles[odir]); \
    loop_free(last); \
    current = last = unshuffled; \
  } \
  return last; \
}

GENERIC_EXCHANGE(unsigned, uints)
GENERIC_EXCHANGE(unsigned long, ulongs)
GENERIC_EXCHANGE(double, doubles)

void free_exchanger(struct exchanger* ex)
{
  for (unsigned i = 0; i < 2; ++i) {
    comm_free(ex->comms[i]);
    loop_free(ex->ranks[i]);
    loop_free(ex->msg_counts[i]);
    loop_free(ex->msg_offsets[i]);
    loop_free(ex->shuffles[i]);
    loop_free(ex->msg_of_items[i]);
    loop_free(ex->items_of_roots_offsets[i]);
  }
  loop_host_free(ex);
}

void exchange_tag(struct exchanger* ex, struct const_tag* t,
    struct tags* into, enum exch_dir dir, enum exch_start start)
{
  void* data_out = 0;
  switch (t->type) {
    case TAG_U8:
      break;
    case TAG_U32:
      data_out = exchange_uints(ex, t->ncomps, t->d.u32, dir, start);
      break;
    case TAG_U64:
      data_out = exchange_ulongs(ex, t->ncomps, t->d.u64, dir, start);
      break;
    case TAG_F64:
      data_out = exchange_doubles(ex, t->ncomps, t->d.f64, dir, start);
      break;
  }
  if (find_tag(into, t->name))
    modify_tag(into, t->name, data_out);
  else
    add_tag(into, t->type, t->name, t->ncomps, data_out);
}

void exchange_tags(struct exchanger* ex, struct tags* from,
    struct tags* into, enum exch_dir dir, enum exch_start start)
{
  for (unsigned i = 0; i < count_tags(from); ++i)
    exchange_tag(ex, get_tag(from, i), into, dir, start);
}
