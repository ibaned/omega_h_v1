#include "owners_from_verts.h"

#include "comm.h"
#include "exchanger.h"
#include "global.h"
#include "ints.h"
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
    unsigned const* spoke_ranks_of_ecopies,
    unsigned const* spoke_idxs_of_ecopies,
    unsigned e,
    unsigned e2)
{
  for (unsigned i = 0; i < (verts_per_elem - 1); ++i) {
    if (spoke_ranks_of_ecopies[e  * (verts_per_elem - 1) + i] !=
        spoke_ranks_of_ecopies[e2 * (verts_per_elem - 1) + i])
      return 0;
    if (spoke_idxs_of_ecopies[e  * (verts_per_elem - 1) + i] !=
        spoke_idxs_of_ecopies[e2 * (verts_per_elem - 1) + i])
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
  unsigned* hub_rank_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* hub_idx_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* spoke_ranks_of_elems = LOOP_MALLOC(unsigned,
      nelems * (verts_per_elem - 1));
  unsigned* spoke_idxs_of_elems = LOOP_MALLOC(unsigned,
      nelems * (verts_per_elem - 1));
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    unsigned hub_vert = verts_of_elem[0];
    hub_rank_of_elems[i] = own_rank_of_verts[hub_vert];
    hub_idx_of_elems[i] = own_idx_of_verts[hub_vert];
    unsigned* spoke_parts = spoke_ranks_of_elems +
      i * (verts_per_elem - 1);
    unsigned* spoke_idxs = spoke_idxs_of_elems +
      i * (verts_per_elem - 1);
    for (unsigned j = 1; j < verts_per_elem; ++j) {
      spoke_parts[j - 1] = own_rank_of_verts[verts_of_elem[j]];
      spoke_idxs[j - 1] = own_idx_of_verts[verts_of_elem[j]];
    }
  }
  struct exchanger* ex = new_exchanger(nelems, hub_rank_of_elems);
  loop_free(hub_rank_of_elems);
  set_exchanger_dests(ex, nverts, hub_idx_of_elems);
  loop_free(hub_idx_of_elems);
  unsigned* elem_idxs = uints_linear(nelems);
  unsigned* idx_of_ecopies = exchange_uints(ex, 1, elem_idxs, EX_FOR, EX_ITEM);
  loop_free(elem_idxs);
  unsigned* spoke_ranks_of_ecopies = exchange_uints(ex, verts_per_elem - 1,
      spoke_ranks_of_elems, EX_FOR, EX_ITEM);
  loop_free(spoke_ranks_of_elems);
  unsigned* spoke_idxs_of_ecopies = exchange_uints(ex, verts_per_elem - 1,
      spoke_idxs_of_elems, EX_FOR, EX_ITEM);
  loop_free(spoke_idxs_of_elems);
  unsigned* recv_nelems = LOOP_MALLOC(unsigned, ex->nmsgs[EX_REV]);
  comm_sync_uint(ex->comms[EX_FOR], nelems, recv_nelems);
  unsigned necopies = ex->nitems[EX_REV];
  unsigned* own_ranks = LOOP_MALLOC(unsigned, necopies);
  unsigned* own_idxs = LOOP_MALLOC(unsigned, necopies);
  for (unsigned i = 0; i < necopies; ++i)
    own_ranks[i] = INVALID;
  /* okay, this algorithm is a mess.
     we have all *copies* of elements which call
     this vertex their hub vertex, and we simultaneously
     have to keep track of which copies are actually
     copies of the same element and also assign
     an owner to that subset of the copies */
  unsigned const* ecopies_of_verts_offsets =
    ex->items_of_roots_offsets[EX_REV];
  unsigned const* recv_of_ecopies = ex->msg_of_items[EX_REV];
  unsigned const* recv_ranks = ex->ranks[EX_REV];
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first = ecopies_of_verts_offsets[i];
    unsigned end = ecopies_of_verts_offsets[i + 1];
    for (unsigned j = first; j < end; ++j) {
      if (own_ranks[j] != INVALID)
        continue; /* an owner was already chosen */
      unsigned own_recv = recv_of_ecopies[j];
      unsigned own_idx = idx_of_ecopies[j];
      for (unsigned k = j + 1; k < end; ++k) {
        if (!are_same_element(verts_per_elem,
              spoke_ranks_of_ecopies,
              spoke_idxs_of_ecopies,
              j, k))
          continue;
        unsigned recv2 = recv_of_ecopies[k];
        unsigned idx2 = idx_of_ecopies[k];
        if ((recv_nelems[recv2] < recv_nelems[own_recv]) ||
            ((recv_nelems[recv2] == recv_nelems[own_recv]) &&
             (recv_ranks[recv2] < recv_ranks[own_recv]))) {
          own_recv = recv2;
          own_idx = idx2;
        }
      }
      for (unsigned k = j + 1; k < end; ++k) {
        if (!are_same_element(verts_per_elem,
              spoke_ranks_of_ecopies,
              spoke_idxs_of_ecopies,
              j, k))
          continue;
        own_ranks[k] = recv_ranks[own_recv];
        own_idxs[k] = own_idx;
      }
    }
  }
  loop_free(recv_nelems);
  loop_free(spoke_ranks_of_ecopies);
  loop_free(spoke_idxs_of_ecopies);
  loop_free(idx_of_ecopies);
  *p_own_ranks = exchange_uints(ex, 1, own_ranks, EX_REV, EX_ITEM);
  *p_own_idxs = exchange_uints(ex, 1, own_idxs, EX_REV, EX_ITEM);
  loop_free(own_ranks);
  loop_free(own_idxs);
  free_exchanger(ex);
}
