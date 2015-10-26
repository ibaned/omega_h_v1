
#include "comm.h"
#include "global.h"
#include "ints.h"
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
  unsigned* sendbuf = sort_uints_by_category(
      lin_idxs, 1, n, sent_sends, sent_idxs, send_offsets);
  unsigned* recvbuf = LOOP_MALLOC(unsigned, nrecvd);
  comm_exch_uints(gc, 1, sendbuf, send_counts, send_offsets,
      recvbuf, recv_counts, recv_offsets);
  comm_free(gc);
}
