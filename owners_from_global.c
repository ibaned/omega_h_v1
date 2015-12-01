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
  struct exchanger* ex = new_exchanger(n, linsize, lin_ranks, lin_idxs);
  loop_free(lin_ranks);
  loop_free(lin_idxs);
  unsigned const* recv_of_recvd = ex->recv_of_recvd;
  unsigned const* recvd_of_lins = ex->recvd_of_dests;
  unsigned const* recvd_of_lin_offsets = ex->recvd_of_dests_offsets;
  unsigned* orig_idxs = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    orig_idxs[i] = i;
  unsigned* orig_idxs_recvd = exchange_uints(ex, 1, orig_idxs);
  loop_free(orig_idxs);
  unsigned* own_rank_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  unsigned* own_idx_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  unsigned* recv_nents = LOOP_MALLOC(unsigned, ex->nrecvs);
  comm_sync_uint(ex->forward_comm, n, recv_nents);
  for (unsigned i = 0; i < linsize; ++i) {
    unsigned first = recvd_of_lin_offsets[i];
    unsigned end = recvd_of_lin_offsets[i + 1];
    unsigned recvd = recvd_of_lins[first];
    unsigned own_recv = recv_of_recvd[recvd];
    unsigned own_idx = orig_idxs_recvd[recvd];
    for (unsigned j = first + 1; j < end; ++j) {
      recvd = recvd_of_lins[j];
      unsigned recv2 = ex->recv_ranks[recv_of_recvd[recvd]];
      unsigned idx2 = orig_idxs_recvd[recvd];
      if ((recv_nents[recv2] < recv_nents[own_recv]) ||
          ((recv_nents[recv2] == recv_nents[own_recv]) &&
           (ex->recv_ranks[recv2] < ex->recv_ranks[own_recv]))) {
        own_recv = recv2;
        own_idx = idx2;
      }
    }
    for (unsigned j = first; j < end; ++j) {
      recvd = recvd_of_lins[j];
      own_rank_of_recvd[recvd] = ex->recv_ranks[own_recv];
      own_idx_of_recvd[recvd] = own_idx;
    }
  }
  loop_free(orig_idxs_recvd);
  loop_free(recv_nents);
  *p_own_ranks = unexchange_uints(ex, 1, own_rank_of_recvd);
  *p_own_idxs = unexchange_uints(ex, 1, own_idx_of_recvd);
  loop_free(own_rank_of_recvd);
  loop_free(own_idx_of_recvd);
  free_exchanger(ex);
}
