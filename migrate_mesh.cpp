#include "migrate_mesh.hpp"

#include <cassert>

#include "close_partition.hpp"
#include "exchanger.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "tables.hpp"

namespace omega_h {

LOOP_KERNEL(use_lids_kern,
    unsigned const* use_offsets,
    unsigned const* copy_offsets,
    unsigned const* use_msgs,
    unsigned const* use_msg_ranks,
    unsigned const* copy_msgs,
    unsigned const* copy_msg_ranks,
    unsigned const* copy_lids,
    unsigned* use_lids)
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
  LOOP_EXEC(use_lids_kern, nowners,
      use_offsets,
      copy_offsets,
      use_msgs,
      use_msg_ranks,
      copy_msgs,
      copy_msg_ranks,
      copy_lids,
      use_lids);
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
  unsigned* verts_of_ents = exchange(use_to_own, 1, lid_of_uses,
      EX_REV, EX_ITEM);
  loop_free(lid_of_uses);
  return verts_of_ents;
}

static struct exchanger* close_uses_from_above(
    struct mesh* m,
    unsigned high_dim,
    unsigned low_dim,
    struct exchanger* high_push)
{
  unsigned* use_own_rank_sent;
  unsigned* use_own_id_sent;
  unsigned* high_use_offsets_sent;
  get_down_use_owners(m, high_dim, low_dim,
      &use_own_rank_sent, &use_own_id_sent, &high_use_offsets_sent);
  unsigned* use_own_rank_recvd;
  unsigned* use_own_id_recvd;
  unsigned* high_use_offsets_recvd;
  push_use_owners(high_push,
      use_own_rank_sent, use_own_id_sent, high_use_offsets_sent,
      &use_own_rank_recvd, &use_own_id_recvd, &high_use_offsets_recvd);
  loop_free(use_own_rank_sent);
  loop_free(use_own_id_sent);
  loop_free(high_use_offsets_sent);
  loop_free(high_use_offsets_recvd);
  unsigned nhighs_recvd = high_push->nitems[EX_REV];
  unsigned nlows = mesh_count(m, low_dim);
  unsigned nlows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned nuses = nhighs_recvd * nlows_per_high;
  struct exchanger* use_to_own = new_exchanger(nuses, use_own_rank_recvd);
  loop_free(use_own_rank_recvd);
  set_exchanger_dests(use_to_own, nlows, use_own_id_recvd);
  loop_free(use_own_id_recvd);
  return use_to_own;
}

static void move_ents(struct mesh* m, struct mesh* m_out, unsigned dim,
    unsigned nents, unsigned* verts_of_ents,
    struct exchanger* push)
{
  mesh_set_ents(m_out, dim, nents, verts_of_ents);
  mesh_tag_globals(m, dim);
  push_tags(push, mesh_tags(m, dim), mesh_tags(m_out, dim));
  mesh_parallel_from_tags(m_out, dim);
}

void migrate_mesh(struct mesh* m, struct exchanger* elem_push)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems_recvd = elem_push->nitems[EX_REV];
  struct exchanger* vert_use_to_own = close_uses_from_above(
      m, dim, 0, elem_push);
  struct exchanger* vert_push = close_partition_exchanger(vert_use_to_own);
  unsigned nverts_recvd = vert_push->nitems[EX_REV];
  struct mesh* m_out = new_mesh(dim, mesh_get_rep(m), 1);
  move_ents(m, m_out, 0, nverts_recvd, 0, vert_push);
  unsigned* lids = uints_linear(nverts_recvd, 1);
  unsigned* lid_of_copies = exchange(vert_push, 1, lids,
      EX_REV, EX_ITEM);
  loop_free(lids);
  unsigned* verts_of_elems = push_connectivity(
      vert_use_to_own, vert_push, lid_of_copies);
  free_exchanger(vert_use_to_own);
  move_ents(m, m_out, dim, nelems_recvd, verts_of_elems, elem_push);
  for (unsigned d = 1; d < dim; ++d) {
    if (!mesh_has_dim(m, d))
      continue;
    struct exchanger* ent_use_to_own = close_uses_from_above(
        m, dim, d, elem_push);
    struct exchanger* ent_push = close_partition_exchanger(ent_use_to_own);
    free_exchanger(ent_use_to_own);
    vert_use_to_own = close_uses_from_above(
        m, d, 0, ent_push);
    unsigned* verts_of_ents = push_connectivity(
        vert_use_to_own, vert_push, lid_of_copies);
    free_exchanger(vert_use_to_own);
    unsigned nents_recvd = ent_push->nitems[EX_REV];
    move_ents(m, m_out, d, nents_recvd, verts_of_ents, ent_push);
    free_exchanger(ent_push);
  }
  loop_free(lid_of_copies);
  free_exchanger(vert_push);
  overwrite_mesh(m, m_out);
}

void migrate_mesh(
    struct mesh* m,
    unsigned nelems_recvd,
    unsigned const* recvd_elem_ranks,
    unsigned const* recvd_elem_ids)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  struct exchanger* elem_push = make_reverse_exchanger(nelems,
      nelems_recvd, recvd_elem_ranks, recvd_elem_ids);
  migrate_mesh(m, elem_push);
  free_exchanger(elem_push);
}

}
