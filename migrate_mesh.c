#include "migrate_mesh.h"

#include <assert.h>

#include "close_partition.h"
#include "exchanger.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"

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

static unsigned* push_connectivity(
    /* exhanger from vertex uses of entities in new mesh
       to vertex owners in the old mesh */
    struct exchanger* use_to_own,
    /* exchanger from entity owners in old mesh to entity copies in new mesh */
    struct exchanger* ent_push,
    /* local IDs in the new mesh of copies of vertices, organized
       by their old mesh owners on this part */
    unsigned const* lid_of_copies)
{
  unsigned* lid_of_uses = use_lids_from_copy_lids(
      use_to_own, ent_push, lid_of_copies);
  unsigned* verts_of_ents = exchange_uints(use_to_own, 1, lid_of_uses,
      EX_REV, EX_ITEM);
  loop_free(lid_of_uses);
  return verts_of_ents;
}

void migrate_mesh(
    struct mesh** p_m,
    unsigned nelems_recvd,
    unsigned const* recvd_elem_ranks,
    unsigned const* recvd_elem_ids)
{
  struct mesh* m = *p_m;
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  unsigned nverts_per_elem = the_down_degrees[dim][0];
  unsigned* use_own_rank_sent;
  unsigned* use_own_id_sent;
  unsigned* elem_use_offsets_sent;
  get_down_use_owners(m, dim, 0, &use_own_rank_sent, &use_own_id_sent,
      &elem_use_offsets_sent);
  struct exchanger* elem_pull = new_exchanger(nelems_recvd, recvd_elem_ranks);
  set_exchanger_dests(elem_pull, nelems, recvd_elem_ids);
  unsigned* use_own_rank_recvd;
  unsigned* use_own_id_recvd;
  unsigned* elem_use_offsets_recvd;
  pull_use_owners(elem_pull,
      use_own_rank_sent, use_own_id_sent, elem_use_offsets_sent,
      &use_own_rank_recvd, &use_own_id_recvd, &elem_use_offsets_recvd);
  loop_free(use_own_rank_sent);
  loop_free(use_own_id_sent);
  loop_free(elem_use_offsets_sent);
  loop_free(elem_use_offsets_recvd);
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
  unsigned* verts_of_elems = push_connectivity(
      use_to_own, vert_push, lid_of_copies);
  loop_free(lid_of_copies);
  free_exchanger(use_to_own);
  struct mesh* m_out = new_mesh(dim);
  mesh_set_ents(m_out, 0, nverts_recvd, 0);
  mesh_set_ents(m_out, dim, nelems_recvd, verts_of_elems);
  exchange_tags(vert_push, mesh_tags(m, 0), mesh_tags(m_out, 0),
      EX_FOR, EX_ROOT);
  free_exchanger(vert_push);
  exchange_tags(elem_pull, mesh_tags(m, dim), mesh_tags(m_out, dim),
      EX_REV, EX_ROOT);
  free_exchanger(elem_pull);
  free_mesh(*p_m);
  *p_m = m_out;
}
