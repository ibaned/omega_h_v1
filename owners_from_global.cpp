#include "owners_from_global.hpp"

#include <cassert>

#include "comm.hpp"
#include "exchanger.hpp"
#include "global.hpp"
#include "ints.hpp"
#include "invert_map.hpp"
#include "loop.hpp"
#include "tables.hpp"

namespace omega_h {

static void setup_linpart(
    unsigned n,
    unsigned long const* global_in,
    struct exchanger** p_ex,
    unsigned** p_orig_idxs_recvd)
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
  unsigned* orig_idxs = uints_linear(n, 1);
  unsigned* orig_idxs_recvd = exchange(ex, 1, orig_idxs, EX_FOR, EX_ITEM);
  loop_free(orig_idxs);
  *p_ex = ex;
  *p_orig_idxs_recvd = orig_idxs_recvd;
}

void owners_from_global(
    unsigned n,
    unsigned long const* global_in,
    unsigned** p_own_ranks,
    unsigned** p_own_idxs)
{
  struct exchanger* ex;
  unsigned* orig_idxs_recvd;
  setup_linpart(n, global_in, &ex, &orig_idxs_recvd);
  unsigned linsize = ex->nroots[EX_REV];
  unsigned const* recvd_of_lin_offsets =
    ex->items_of_roots_offsets[EX_REV];
  unsigned const* recv_of_recvd = ex->msg_of_items[EX_REV];
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
      unsigned recv2 = recv_of_recvd[j];
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
  *p_own_ranks = exchange(ex, 1, own_rank_of_recvd, EX_REV, EX_ITEM);
  *p_own_idxs = exchange(ex, 1, own_idx_of_recvd, EX_REV, EX_ITEM);
  loop_free(own_rank_of_recvd);
  loop_free(own_idx_of_recvd);
  free_exchanger(ex);
}

void own_idxs_from_global(
    unsigned n,
    unsigned long const* global_in,
    unsigned const* own_ranks_in,
    unsigned** p_own_idxs)
{
  struct exchanger* ex;
  unsigned* orig_idxs_recvd;
  setup_linpart(n, global_in, &ex, &orig_idxs_recvd);
  unsigned linsize = ex->nroots[EX_REV];
  unsigned const* recvd_of_lin_offsets =
    ex->items_of_roots_offsets[EX_REV];
  unsigned const* recv_of_recvd = ex->msg_of_items[EX_REV];
  unsigned* own_rank_of_recvd = exchange(ex, 1, own_ranks_in,
      EX_FOR, EX_ITEM);
  unsigned nrecvd = ex->nitems[EX_REV];
  unsigned* own_idx_of_recvd = LOOP_MALLOC(unsigned, nrecvd);
  unsigned const* recv_ranks = ex->ranks[EX_REV];
  for (unsigned i = 0; i < linsize; ++i) {
    unsigned first = recvd_of_lin_offsets[i];
    unsigned end = recvd_of_lin_offsets[i + 1];
    unsigned own_rank = own_rank_of_recvd[first];
    unsigned own_idx = INVALID;
    for (unsigned j = first; j < end; ++j) {
      assert(own_rank_of_recvd[j] == own_rank);
      if (recv_ranks[recv_of_recvd[j]] == own_rank) {
        own_idx = orig_idxs_recvd[j];
        break;
      }
    }
    assert(own_idx != INVALID);
    for (unsigned j = first; j < end; ++j) {
      own_idx_of_recvd[j] = own_idx;
    }
  }
  loop_free(own_rank_of_recvd);
  loop_free(orig_idxs_recvd);
  *p_own_idxs = exchange(ex, 1, own_idx_of_recvd, EX_REV, EX_ITEM);
  loop_free(own_idx_of_recvd);
  free_exchanger(ex);
}

}
