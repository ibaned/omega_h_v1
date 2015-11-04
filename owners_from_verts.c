#include "owners_from_verts.h"

#include "global.h"
#include "loop.h"
#include "tables.h"

void owners_from_verts(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    unsigned const* own_part_of_verts,
    unsigned const* own_idx_of_verts,
    unsigned** p_own_parts,
    unsigned** p_own_idxs)
{
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned nuses = nelems * verts_per_elem;
  unsigned* own_vert_part_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* own_vert_idx_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* vert_own_parts_of_elems = LOOP_MALLOC(unsigned,
      nelems * (verts_per_elem - 1));
  unsigned* vert_own_idxs_of_elems = LOOP_MALLOC(unsigned,
      nelems * (verts_per_elem - 1));
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    unsigned own_vert = verts_of_elem[0];
    unsigned own_vert_part_of_elems[i] = own_part_of_verts[own_vert];
    unsigned own_vert_idx_of_elems[i] = own_idx_of_verts[own_vert];
    unsigned* vert_own_parts = vert_own_parts_of_elems +
      i * (verts_per_elem - 1);
    unsigned* vert_own_idxs = vert_own_idxs_of_elems +
      i * (verts_per_elem - 1);
    for (unsigned j = 1; j < verts_per_elem; ++j) {
      vert_own_parts[j - 1] = own_part_of_verts[verts_of_elem[j]];
      vert_own_idxs[j - 1] = own_idx_of_verts[verts_of_elem[j]];
    }
  }
  /* communication setup */
  unsigned* elem_sends;
  unsigned* elem_send_idxs;
  unsigned nsends;
  unsigned* send_parts;
  unsigned* send_counts;
  categorize_by_part(own_vert_part_of_elems, nelems,
      &elem_sends, &elem_send_idxs,
      &nsends, &send_parts, &send_counts);
  unsigned* send_offsets = uints_exscan(send_counts, nsends);
  struct comm* c = comm_graph(comm_using(), nsends, send_parts, send_counts);
  unsigned nrecvs;
  unsigned* recv_parts;
  unsigned* recv_counts;
  comm_recvs(c, &nrecvs, &recv_parts, &recv_counts);
  unsigned* recv_offsets = uints_exscan(recv_counts, nrecvs);
  unsigned* orig_idxs = LOOP_MALLOC(unsigned, nelems);
  for (unsigned i = 0; i < nelems; ++i)
    orig_idxs[i] = i;
  /* send original indices */
  unsigned* sent_orig_idxs = sort_uints_by_category(
      orig_idxs, 1, nelems, elem_sends, elem_send_idxs, send_offsets);
  unsigned* recvd_orig_idxs = comm_exch_uints(c, 1,
      sent_orig_idxs, send_counts, send_offsets, recv_counts, recv_offsets);
  /* send own_vert_part_of_elems */
  unsigned* sent_ovpe = sort_uints_by_category(
      own_vert_part_of_elems, 1, nelems, elem_sends, elem_send_idxs, send_offsets);
  unsigned* recvd_ovpe = comm_exch_uints(c, 1,
      sent_ovpe, send_counts, send_offsets, recv_counts, recv_offsets);
  /* send own_vert_idx_of_elems */
  unsigned* sent_ovie = sort_uints_by_category(
      own_vert_idx_of_elems, 1, nelems, elem_sends, elem_send_idxs, send_offsets);
  unsigned* recvd_ovie = comm_exch_uints(c, 1,
      sent_ovie, send_counts, send_offsets, recv_counts, recv_offsets);
}
