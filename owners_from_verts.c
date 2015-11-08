#include "owners_from_verts.h"

#include "exchanger.h"
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
  struct exchanger* ex = new_exchanger(nelems, nverts,
      own_vert_part_of_elems, own_vert_idx_of_elems);
  unsigned* orig_idxs = LOOP_MALLOC(unsigned, nelems);
  for (unsigned i = 0; i < nelems; ++i)
    orig_idxs[i] = i;
  unsigned* recvd_orig_idxs = exchange_uints(ex, 1, orig_idxs);
  unsigned* vop_of_recvd = exchange_uints(ex, verts_per_elem - 1,
      vert_own_parts_of_elems);
  unsigned* voi_of_recvd = exchange_uints(ex, verts_per_elem - 1,
      vert_own_idx_of_elems);
  unsigned* recv_nelems = LOOP_MALLOC(unsigned, ex->nrecvs);
  comm_sync_uint(ex->forward_comm, nelems, recv_nelems);
  unsigned* op_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  unsigned* oi_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  for (unsigned i = 0; i < ex->nrecvd; ++i)
    op_of_recvd[i] = INVALID;
  /* okay, this algorithm is a mess.
     we have all *copies* of elements which call
     this vertex their vertex #0, and we simultaneously
     have to keep track of which copies are actually
     copies of the same element and also assign
     an owner to that subset of the copies */
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first = ex->recvd_of_dests_offsets[i];
    unsigned end = ex->recvd_of_dests_offsets[i + 1];
    for (unsigned j = first; j < end; ++j) {
      unsigned elem = ex->recvd_of_dests[j];
      if (op_of_recvd[elem] != INVALID)
        continue; /* an owner was already chosen */
      unsigned op = ex->recv_ranks
      for (unsigned k = first; k < end; ++k) {
      }
    }
  }
}
