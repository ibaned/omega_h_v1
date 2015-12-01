#include "owners_from_global.h"

#include "comm.h"
#include "exchanger.h"
#include "global.h"
#include "ints.h"
#include "invert_map.h"
#include "loop.h"

void owners_from_global(
    unsigned n,
    unsigned long const* global_in,
    unsigned** p_own_ranks,
    unsigned** p_own_idxs)
{
  unsigned long total = comm_max_ulong(ulongs_max(global_in, n)) + 1;
  unsigned nranks = comm_size();
  unsigned* lin_ranks;
  unsigned* lin_idxs;
  global_to_linpart(global_in, n, total, nranks, &lin_ranks, &lin_idxs);
  unsigned linsize = linpart_size(total, nranks, comm_rank());
  struct exchanger* ex = new_exchanger(n, lin_ranks);
  loop_free(lin_ranks);
  set_exchanger_dests(ex, linsize, lin_idxs);
  loop_free(lin_idxs);
  unsigned const* recv_of_recvd = ex->msg_of_items[EX_REV];
  unsigned const* recvd_of_lin_offsets =
    ex->items_of_roots_offsets[EX_REV];
  unsigned* orig_idxs = uints_linear(n);
  unsigned* orig_idxs_recvd = exchange_uints(ex, 1, orig_idxs, EX_FOR, EX_ITEM);
  loop_free(orig_idxs);
  unsigned nrecvd = ex->nitems[EX_REV];
  unsigned nrecvs = ex->nmsgs[EX_REV];
  unsigned* own_rank_of_recvd = LOOP_MALLOC(unsigned, nrecvd);
  unsigned* own_idx_of_recvd = LOOP_MALLOC(unsigned, nrecvd);
  unsigned* recv_nents = LOOP_MALLOC(unsigned, nrecvs);
  comm_sync_uint(ex->comms[EX_FOR], n, recv_nents);
  unsigned const* recv_ranks = ex->ranks[EX_REV];
  for (unsigned i = 0; i < linsize; ++i) {
    unsigned first = recvd_of_lin_offsets[i];
    unsigned end = recvd_of_lin_offsets[i + 1];
    unsigned own_recv = recv_of_recvd[first];
    unsigned own_idx = orig_idxs_recvd[first];
    for (unsigned j = first + 1; j < end; ++j) {
      unsigned recv2 = recv_ranks[recv_of_recvd[j]];
      unsigned idx2 = orig_idxs_recvd[j];
      if ((recv_nents[recv2] < recv_nents[own_recv]) ||
          ((recv_nents[recv2] == recv_nents[own_recv]) &&
           (recv_ranks[recv2] < recv_ranks[own_recv]))) {
        own_recv = recv2;
        own_idx = idx2;
      }
    }
    for (unsigned j = first; j < end; ++j) {
      own_rank_of_recvd[j] = recv_ranks[own_recv];
      own_idx_of_recvd[j] = own_idx;
    }
  }
  loop_free(orig_idxs_recvd);
  loop_free(recv_nents);
  *p_own_ranks = exchange_uints(ex, 1, own_rank_of_recvd, EX_REV, EX_ITEM);
  *p_own_idxs = exchange_uints(ex, 1, own_idx_of_recvd, EX_REV, EX_ITEM);
  loop_free(own_rank_of_recvd);
  loop_free(own_idx_of_recvd);
  free_exchanger(ex);
}
