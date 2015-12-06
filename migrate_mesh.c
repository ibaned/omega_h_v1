#include "migrate_mesh.h"

#include <assert.h>

#include "close_partition.h"
#include "exchanger.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"

/* for each element in the mesh, for each
   vertex it uses, return the owner of that
   vertex. the result is of size (nelems * verts_per_elem) */

static void get_vert_use_owners_of_elems(
    struct mesh* m,
    /* own rank of vertex uses of elements */
    unsigned** p_use_own_ranks,
    unsigned** p_use_own_ids)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  unsigned nverts_per_elem = the_down_degrees[dim][0];
  unsigned const* verts_of_elems = mesh_ask_down(m, dim, 0);
  unsigned const* vert_own_ranks = mesh_ask_own_ranks(m, 0);
  unsigned const* vert_own_ids = mesh_ask_own_ids(m, 0);
  unsigned* use_own_ranks = LOOP_MALLOC(unsigned,
      nelems * nverts_per_elem);
  unsigned* use_own_ids = LOOP_MALLOC(unsigned,
      nelems * nverts_per_elem);
  for (unsigned i = 0; i < nelems; ++i) {
    for (unsigned j = 0; j < nverts_per_elem; ++j) {
      unsigned vert = verts_of_elems[i * nverts_per_elem + j];
      use_own_ranks[i * nverts_per_elem + j] = vert_own_ranks[vert];
      use_own_ids[i * nverts_per_elem + j] = vert_own_ids[vert];
    }
  }
  *p_use_own_ranks = use_own_ranks;
  *p_use_own_ids = use_own_ids;
}

static unsigned* use_lids_from_copy_lids(
    struct exchanger* use_to_own,
    struct exchanger* own_to_copy,
    unsigned const* copy_lids)
{
  unsigned nuses = use_to_own->nitems[EX_REV];
  unsigned nowners = own_to_copy->nroots[EX_FOR];
  unsigned const* use_offsets =
    use_to_own->items_of_roots_offsets[EX_REV];
  unsigned const* use_msgs =
    use_to_own->msg_of_items[EX_REV];
  unsigned const* use_msg_ranks =
    use_to_own->ranks[EX_REV];
  unsigned const* copy_offsets =
    own_to_copy->items_of_roots_offsets[EX_FOR];
  unsigned const* copy_msgs =
    own_to_copy->msg_of_items[EX_FOR];
  unsigned const* copy_msg_ranks =
    own_to_copy->ranks[EX_FOR];
  unsigned* use_lids = LOOP_MALLOC(unsigned, nuses);
  for (unsigned i = 0; i < nowners; ++i) {
    unsigned fu = use_offsets[i];
    unsigned eu = use_offsets[i + 1];
    unsigned fc = copy_offsets[i];
    unsigned ec = copy_offsets[i + 1];
    for (unsigned j = fu; j < eu; ++j) {
      unsigned msg = use_msgs[j];
      unsigned rank = use_msg_ranks[msg];
      unsigned lid = INVALID;
      for (unsigned k = fc; k < ec; ++k)
        if (rank == copy_msg_ranks[copy_msgs[k]]) {
          lid = copy_lids[k];
          break;
        }
      assert(lid != INVALID);
      use_lids[j] = lid;
    }
  }
  return use_lids;
}

static struct mesh* migrate_element_topology(
    struct mesh* m,
    unsigned nelems_recvd,
    unsigned const* recvd_elem_ranks,
    unsigned const* recvd_elem_ids,
    /* this exchanger sends from new elements to old elements */
    struct exchanger** p_elem_pull,
    /* this exchanger sends from old vertex owners to new vertex copies */
    struct exchanger** p_vert_push)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  unsigned nverts_per_elem = the_down_degrees[dim][0];
  unsigned* use_own_rank_sent;
  unsigned* use_own_id_sent;
  get_vert_use_owners_of_elems(m, &use_own_rank_sent, &use_own_id_sent);
  struct exchanger* elem_pull = new_exchanger(nelems_recvd, recvd_elem_ranks);
  set_exchanger_dests(elem_pull, nelems, recvd_elem_ids);
  /* own rank of vertex uses of received elements */
  unsigned* use_own_rank_recvd = exchange_uints(elem_pull, nverts_per_elem,
      use_own_rank_sent, EX_REV, EX_ROOT);
  loop_free(use_own_rank_sent);
  unsigned* use_own_id_recvd = exchange_uints(elem_pull, nverts_per_elem,
      use_own_id_sent, EX_REV, EX_ROOT);
  loop_free(use_own_id_sent);
  unsigned nverts = mesh_count(m, 0);
  struct exchanger* use_to_own;
  struct exchanger* vert_push;
  unsigned* uses_by_copies_offsets = uints_linear(nelems_recvd, nverts_per_elem);
  close_partition_exchangers(nelems_recvd, nverts, uses_by_copies_offsets,
      use_own_rank_recvd, use_own_id_recvd, &use_to_own, &vert_push);
  loop_free(use_own_rank_recvd);
  loop_free(use_own_id_recvd);
  loop_free(uses_by_copies_offsets);
  unsigned nverts_recvd = vert_push->nitems[EX_REV];
  unsigned* lids = uints_linear(nverts_recvd, 1);
  unsigned* lid_of_copies = exchange_uints(vert_push, 1, lids,
      EX_REV, EX_ITEM);
  loop_free(lids);
  unsigned* lid_of_rvus = use_lids_from_copy_lids(
      use_to_own, vert_push, lid_of_copies);
  loop_free(lid_of_copies);
  unsigned* verts_of_elems = exchange_uints(use_to_own, 1, lid_of_rvus,
      EX_REV, EX_ITEM);
  loop_free(lid_of_rvus);
  free_exchanger(use_to_own);
  struct mesh* m_out = new_mesh(dim);
  mesh_set_ents(m_out, 0, nverts_recvd, 0);
  mesh_set_ents(m_out, dim, nelems_recvd, verts_of_elems);
  *p_elem_pull = elem_pull;
  *p_vert_push = vert_push;
  return m_out;
}

void migrate_mesh(struct mesh** p_m,
    unsigned nelems_recvd,
    unsigned const* recvd_elem_ranks,
    unsigned const* recvd_elem_ids)
{
  struct exchanger* elem_pull;
  struct exchanger* vert_push;
  struct mesh* m_in = *p_m;
  struct mesh* m_out = migrate_element_topology(m_in,
      nelems_recvd, recvd_elem_ranks, recvd_elem_ids,
      &elem_pull, &vert_push);
  exchange_tags(vert_push, mesh_tags(m_in, 0), mesh_tags(m_out, 0),
      EX_FOR, EX_ROOT);
  free_exchanger(vert_push);
  unsigned dim = mesh_dim(m_in);
  exchange_tags(elem_pull, mesh_tags(m_in, dim), mesh_tags(m_out, dim),
      EX_REV, EX_ROOT);
  free_exchanger(elem_pull);
  free_mesh(*p_m);
  *p_m = m_out;
}
