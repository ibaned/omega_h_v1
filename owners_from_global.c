#include "owners_from_global.h"

#include "comm.h"
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
  unsigned long total = comm_add_ulong(ulongs_max(global_in, n));
  unsigned nparts = comm_size();
  unsigned* lin_parts;
  unsigned* lin_idxs;
  global_to_linpart(global_in, n, total, nparts, &lin_parts, &lin_idxs);
  unsigned* sent_sends;
  unsigned* sent_idxs;
  unsigned nsends;
  unsigned* send_parts;
  unsigned* send_counts;
  categorize_by_part(lin_parts, n, &sent_sends, &sent_idxs,
      &nsends, &send_parts, &send_counts);
  struct comm* gc = comm_graph(comm_using(), nsends, send_parts, send_counts);
  unsigned nrecvs;
  unsigned** recv_parts;
  unsigned** recv_counts;
  comm_recvs(gc, &nrecvs, &recv_parts, &recv_counts);
  unsigned* send_offsets = uints_exscan(send_counts, nsends);
  unsigned* recv_offsets = uints_exscan(recv_counts, nrecvs);
  unsigned nrecvd = recv_offsets[nrecvs];
  unsigned* lin_idxs_sent = sort_uints_by_category(
      lin_idxs, 1, n, sent_sends, sent_idxs, send_offsets);
  unsigned* lin_idxs_recvd = LOOP_MALLOC(unsigned, nrecvd);
  comm_exch_uints(gc, 1, lin_idxs_sent, send_counts, send_offsets,
      lin_idxs_recvd, recv_counts, recv_offsets);
  unsigned* recv_of_recvd = LOOP_MALLOC(unsigned, nrecvd);
  for (unsigned i = 0; i < nrecvs; ++i) {
    unsigned first = recv_offsets[i];
    unsigned end = recv_offsets[i + 1];
    for (unsigned j = first; j < end; ++j)
      recv_of_recvd[j] = i;
  }
  unsigned linsize = linpart_size(total, nparts, comm_rank());
  unsigned* recvd_of_lins;
  unsigned* recvd_of_lin_offsets;
  invert_map(nrecvd, lin_idxs_recvd, linsize,
      &recvd_of_lins, &recvd_of_lin_offsets);
  unsigned* orig_idxs = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    orig_idxs[i] = i;
  unsigned* orig_idxs_sent = sort_uints_by_category(
      orig_idxs, 1, n, sent_sends, sent_idxs, send_offsets);
  unsigned* orig_idxs_recvd = LOOP_MALLOC(unsigned, nrecvd);
  comm_exch_uints(gc, 1, orig_idxs_sent, send_counts, send_offsets,
      orig_idxs_recvd, recv_counts, recv_offsets);
  unsigned* own_part_of_recvd = LOOP_MALLOC(unsigned, nrecvd);
  unsigned* own_idx_of_recvd = LOOP_MALLOC(unsigned, nrecvd);
  unsigned* recv_nents = LOOP_MALLOC(unsigned, nrecvs);
  comm_sync_uint(gc, n, recv_nents);
  comm_free(gc);
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
        own_idx = orig_idxs_recvd[recvd]
      } else if ((recv_nents[recv_of_recvd[recvd]] ==
                  recv_nents[own_recv]) &&
                 (recv_parts[recv_of_recvd[recvd]] <
                  recv_parts[own_recv])) {
        own_recv = recv_of_recvd[recvd];
        own_idx = orig_idxs_recvd[recvd]
      }
    }
    for (unsigned j = first; j < end; ++j) {
      recvd = recvd_of_lins[j];
      own_part_of_recvd[recvd] = recv_parts[own_recv];
      own_idx_of_recvd[recvd] = own_idx;
    }
  }
  struct comm* gci = comm_graph(comm_using(), nrecvs, recv_parts, recv_counts);
  unsigned* own_part_of_sent = LOOP_MALLOC(unsigned, n);
  comm_exch_uints(gci, 1, own_part_of_recvd, recv_counts, recv_offsets,
      own_part_of_sent, send_counts, send_offsets);
  unsigned* own_idx_of_sent = LOOP_MALLOC(unsigned, n);
  comm_exch_uints(gci, 1, own_idx_of_recvd, recv_counts, recv_offsets,
      own_idx_of_sent, send_counts, send_offsets);
  comm_free(gci);
  *p_own_parts = unsort_uints_by_category(
      own_part_of_sent, 1, n, sent_sends, sent_idxs, send_offsets);
  *p_own_idxs = unsort_uints_by_category(
      own_idx_of_sent, 1, n, sent_sends, sent_idxs, send_offsets);
}
