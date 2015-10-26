
#include "comm.h"
#include "global.h"
#include "ints.h"
#include "invert_map.h"
#include "loop.h"

unsigned long* global_renumber(
    unsigned long* global_in,
    unsigned n)
{
  unsigned long total = comm_add_ulong(ulongs_max(global_in, n));
  unsigned nparts = comm_size();
  unsigned* lin_parts;
  unsigned* lin_idxs;
  global_to_linpart(global_in, n, total, nparts, lin_parts, lin_idxs);
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
  unsigned* src_parts_recvd = LOOP_MALLOC(unsigned, nrecvd);
  for (unsigned i = 0; i < nrecvs; ++i) {
    unsigned first = recv_offsets[i];
    unsigned end = recv_offsets[i + 1];
    for (unsigned j = first; j < end; ++j)
      src_parts_recvd[j] = recv_parts[i];
  }
  unsigned linsize = linpart_size(total, nparts, comm_rank());
  unsigned* recvd_of_lins;
  unsigned* recvd_of_lin_offsets;
  invert_map(nrecvd, lin_idxs_recvd, linsize, &recvd_of_lins,
      &recvd_of_lin_offsets);
  comm_free(gc);
}
