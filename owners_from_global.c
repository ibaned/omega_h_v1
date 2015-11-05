#include "owners_from_global.h"

#include "comm.h"
#include "exchanger.h"
#include "global.h"
#include "ints.h"
#include "invert_map.h"
#include "loop.h"

void owners_from_global(
    unsigned n,
    unsigned long* global_in,
    unsigned** p_own_parts,
    unsigned** p_own_idxs)
{
  unsigned long total = comm_max_ulong(ulongs_max(global_in, n)) + 1;
  unsigned nparts = comm_size();
  unsigned* lin_parts;
  unsigned* lin_idxs;
  global_to_linpart(global_in, n, total, nparts, &lin_parts, &lin_idxs);
  struct exchanger* ex = new_exchanger(n, lin_parts);
  loop_free(lin_parts);
  unsigned* lin_idxs_recvd = exchange_uints(ex, 1, lin_idxs);
  loop_free(lin_idxs);
  unsigned* recv_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  for (unsigned i = 0; i < ex->nrecvs; ++i) {
    unsigned first = ex->recv_offsets[i];
    unsigned end = ex->recv_offsets[i + 1];
    for (unsigned j = first; j < end; ++j)
      recv_of_recvd[j] = i;
  }
  unsigned linsize = linpart_size(total, nparts, comm_rank());
  unsigned* recvd_of_lins;
  unsigned* recvd_of_lin_offsets;
  invert_map(ex->nrecvd, lin_idxs_recvd, linsize,
      &recvd_of_lins, &recvd_of_lin_offsets);
  loop_free(lin_idxs_recvd);
  unsigned* orig_idxs = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    orig_idxs[i] = i;
  unsigned* orig_idxs_recvd = exchange_uints(ex, 1, orig_idxs);
  loop_free(orig_idxs);
  unsigned* own_part_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  unsigned* own_idx_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  unsigned* recv_nents = LOOP_MALLOC(unsigned, ex->nrecvs);
  comm_sync_uint(ex->forward_comm, n, recv_nents);
  for (unsigned i = 0; i < linsize; ++i) {
    unsigned first = recvd_of_lin_offsets[i];
    unsigned end = recvd_of_lin_offsets[i + 1];
    unsigned recvd = recvd_of_lins[first];
    unsigned own_recv = recv_of_recvd[first];
    unsigned own_idx = orig_idxs_recvd[first];
    for (unsigned j = first + 1; j < end; ++j) {
      recvd = recvd_of_lins[j];
      if (recv_nents[recv_of_recvd[recvd]] <
          recv_nents[own_recv]) {
        own_recv = recv_of_recvd[recvd];
        own_idx = orig_idxs_recvd[recvd];
      } else if ((recv_nents[recv_of_recvd[recvd]] ==
                  recv_nents[own_recv]) &&
                 (ex->recv_parts[recv_of_recvd[recvd]] <
                  ex->recv_parts[own_recv])) {
        own_recv = recv_of_recvd[recvd];
        own_idx = orig_idxs_recvd[recvd];
      }
    }
    for (unsigned j = first; j < end; ++j) {
      recvd = recvd_of_lins[j];
      own_part_of_recvd[recvd] = ex->recv_parts[own_recv];
      own_idx_of_recvd[recvd] = own_idx;
    }
  }
  loop_free(recv_of_recvd);
  loop_free(recvd_of_lins);
  loop_free(recvd_of_lin_offsets);
  loop_free(orig_idxs_recvd);
  loop_free(recv_nents);
  *p_own_parts = unexchange_uints(ex, 1, own_part_of_recvd);
  *p_own_idxs = unexchange_uints(ex, 1, own_idx_of_recvd);
  loop_free(own_part_of_recvd);
  loop_free(own_idx_of_recvd);
  free_exchanger(ex);
}
