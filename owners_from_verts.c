#include "owners_from_verts.h"

#include "comm.h"
#include "exchanger.h"
#include "global.h"
#include "loop.h"
#include "tables.h"

/* checks whether two received copies (e,e2) describe the same element
   by comparing the owners of their vertices
   (except the first one, which in this case we already
    know is the same)
   we assume (require) that the order of the vertices
   is the same amongst copies */
static unsigned are_same_element(
    unsigned verts_per_elem,
    unsigned const* vert_own_ranks_of_recvd,
    unsigned const* vert_own_idxs_of_recvd,
    unsigned e,
    unsigned e2)
{
  for (unsigned i = 0; i < (verts_per_elem - 1); ++i) {
    if (vert_own_ranks_of_recvd[e  * (verts_per_elem - 1) + i] !=
        vert_own_ranks_of_recvd[e2 * (verts_per_elem - 1) + i])
      return 0;
    if (vert_own_idxs_of_recvd[e  * (verts_per_elem - 1) + i] !=
        vert_own_idxs_of_recvd[e2 * (verts_per_elem - 1) + i])
      return 0;
  }
  return 1;
}

void owners_from_verts(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    unsigned const* own_rank_of_verts,
    unsigned const* own_idx_of_verts,
    unsigned** p_own_ranks,
    unsigned** p_own_idxs)
{
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* own_vert_rank_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* own_vert_idx_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* vert_own_ranks_of_elems = LOOP_MALLOC(unsigned,
      nelems * (verts_per_elem - 1));
  unsigned* vert_own_idxs_of_elems = LOOP_MALLOC(unsigned,
      nelems * (verts_per_elem - 1));
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    unsigned own_vert = verts_of_elem[0];
    own_vert_rank_of_elems[i] = own_rank_of_verts[own_vert];
    own_vert_idx_of_elems[i] = own_idx_of_verts[own_vert];
    unsigned* vert_own_parts = vert_own_ranks_of_elems +
      i * (verts_per_elem - 1);
    unsigned* vert_own_idxs = vert_own_idxs_of_elems +
      i * (verts_per_elem - 1);
    for (unsigned j = 1; j < verts_per_elem; ++j) {
      vert_own_parts[j - 1] = own_rank_of_verts[verts_of_elem[j]];
      vert_own_idxs[j - 1] = own_idx_of_verts[verts_of_elem[j]];
    }
  }
  struct exchanger* ex = new_exchanger(nelems, nverts,
      own_vert_rank_of_elems, own_vert_idx_of_elems);
  loop_free(own_vert_rank_of_elems);
  loop_free(own_vert_idx_of_elems);
  unsigned* orig_idxs = LOOP_MALLOC(unsigned, nelems);
  for (unsigned i = 0; i < nelems; ++i)
    orig_idxs[i] = i;
  unsigned* orig_idx_of_recvd = exchange_uints(ex, 1, orig_idxs);
  loop_free(orig_idxs);
  unsigned* vert_own_ranks_of_recvd = exchange_uints(ex, verts_per_elem - 1,
      vert_own_ranks_of_elems);
  loop_free(vert_own_ranks_of_elems);
  unsigned* vert_own_idxs_of_recvd = exchange_uints(ex, verts_per_elem - 1,
      vert_own_idxs_of_elems);
  loop_free(vert_own_idxs_of_elems);
  unsigned* recv_nelems = LOOP_MALLOC(unsigned, ex->nrecvs);
  comm_sync_uint(ex->forward_comm, nelems, recv_nelems);
  unsigned* own_rank_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  unsigned* own_idx_of_recvd = LOOP_MALLOC(unsigned, ex->nrecvd);
  for (unsigned i = 0; i < ex->nrecvd; ++i)
    own_rank_of_recvd[i] = INVALID;
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
      unsigned e = ex->recvd_of_dests[j];
      if (own_rank_of_recvd[e] != INVALID)
        continue; /* an owner was already chosen */
      unsigned own_recv = ex->recv_of_recvd[e];
      unsigned own_idx = orig_idx_of_recvd[e];
      for (unsigned k = j + 1; k < end; ++k) {
        unsigned e2 = ex->recvd_of_dests[j];
        if (!are_same_element(verts_per_elem,
              vert_own_ranks_of_recvd,
              vert_own_idxs_of_recvd,
              e, e2))
          continue;
        unsigned recv2 = ex->recv_of_recvd[e2];
        unsigned idx2 = orig_idx_of_recvd[e];
        if ((recv_nelems[recv2] < recv_nelems[own_recv]) ||
            ((recv_nelems[recv2] == recv_nelems[own_recv]) &&
             (ex->recv_ranks[recv2] < ex->recv_ranks[own_recv]))) {
          own_recv = recv2;
          own_idx = idx2;
        }
      }
      for (unsigned k = j + 1; k < end; ++k) {
        unsigned e2 = ex->recvd_of_dests[j];
        if (!are_same_element(verts_per_elem,
              vert_own_ranks_of_recvd,
              vert_own_idxs_of_recvd,
              e, e2))
          continue;
        own_rank_of_recvd[e2] = ex->recv_ranks[own_recv];
        own_idx_of_recvd[e2] = own_idx;
      }
    }
  }
  loop_free(recv_nelems);
  loop_free(vert_own_ranks_of_recvd);
  loop_free(vert_own_idxs_of_recvd);
  loop_free(orig_idx_of_recvd);
  /* I think (verts_per_elem - 1) should be 1 here... this code is untested... */
  *p_own_ranks = unexchange_uints(ex, (verts_per_elem - 1), own_rank_of_recvd);
  *p_own_idxs = unexchange_uints(ex, (verts_per_elem - 1), own_idx_of_recvd);
  loop_free(own_idx_of_recvd);
  loop_free(own_rank_of_recvd);
  free_exchanger(ex);
}
