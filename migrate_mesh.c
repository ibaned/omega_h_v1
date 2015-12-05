#include "migrate_mesh.h"

#include <assert.h>

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
    unsigned** p_or_of_vus_of_es,
    unsigned** p_oi_of_vus_of_es)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  unsigned nverts_per_elem = the_down_degrees[dim][0];
  unsigned const* verts_of_elems = mesh_ask_down(m, dim, 0);
  unsigned const* vert_own_ranks = mesh_ask_own_ranks(m, 0);
  unsigned const* vert_own_ids = mesh_ask_own_ids(m, 0);
  unsigned* or_of_vus_of_es = LOOP_MALLOC(unsigned,
      nelems * nverts_per_elem);
  unsigned* oi_of_vus_of_es = LOOP_MALLOC(unsigned,
      nelems * nverts_per_elem);
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * nverts_per_elem;
    unsigned* or_of_vus_of_e = or_of_vus_of_es +
        i * nverts_per_elem;
    unsigned* oi_of_vus_of_e = oi_of_vus_of_es +
        i * nverts_per_elem;
    for (unsigned j = 0; j < nverts_per_elem; ++j) {
      unsigned vert = verts_of_elem[j];
      or_of_vus_of_e[j] = vert_own_ranks[vert];
      oi_of_vus_of_e[j] = vert_own_ids[vert];
    }
  }
  *p_or_of_vus_of_es = or_of_vus_of_es;
  *p_oi_of_vus_of_es = oi_of_vus_of_es;
}

static unsigned get_unique_ranks_of_owner(
    unsigned const* uses_of_owners_offsets,
    unsigned const* msg_of_uses,
    unsigned const* rank_of_msgs,
    unsigned owner,
    unsigned* uranks_of_owner)
{
  unsigned nuranks = 0;
  unsigned first = uses_of_owners_offsets[owner];
  unsigned end = uses_of_owners_offsets[owner + 1];
  for (unsigned j = first; j < end; ++j) {
    unsigned msg = msg_of_uses[j];
    unsigned rank = rank_of_msgs[msg];
    nuranks = add_unique(uranks_of_owner, nuranks, rank);
  }
  return nuranks;
}

static void get_unique_ranks_of_owners(
    unsigned nowners,
    unsigned const* uses_of_owners_offsets,
    unsigned const* msg_of_uses,
    unsigned const* rank_of_msgs,
    unsigned** p_copies_of_owners_offsets,
    unsigned** p_rank_of_copies)
{
  unsigned* ncopies_of_owners = LOOP_MALLOC(unsigned, nowners);
  unsigned nuses = uses_of_owners_offsets[nowners];
  unsigned* scratch = LOOP_MALLOC(unsigned, nuses);
  for (unsigned i = 0; i < nowners; ++i) {
    ncopies_of_owners[i] = get_unique_ranks_of_owner(
        uses_of_owners_offsets, msg_of_uses, rank_of_msgs,
        i, scratch + uses_of_owners_offsets[i]);
  }
  unsigned* copies_of_owners_offsets = uints_exscan(
      ncopies_of_owners, nowners);
  loop_free(ncopies_of_owners);
  unsigned ncopies = copies_of_owners_offsets[nowners];
  unsigned* rank_of_copies = LOOP_MALLOC(unsigned, ncopies);
  for (unsigned i = 0; i < nowners; ++i) {
    unsigned fu = uses_of_owners_offsets[i];
    unsigned fc = copies_of_owners_offsets[i];
    unsigned ec = copies_of_owners_offsets[i + 1];
    unsigned nc = ec - fc;
    for (unsigned j = 0; j < nc; ++j)
      rank_of_copies[fc + j] = scratch[fu + j];
  }
  loop_free(scratch);
  *p_copies_of_owners_offsets = copies_of_owners_offsets;
  *p_rank_of_copies = rank_of_copies;
}

static unsigned* use_lids_from_copy_lids(
    unsigned nuses,
    unsigned nowners,
    unsigned const* uses_of_owners_offsets,
    unsigned const* copies_of_owners_offsets,
    unsigned const* rank_of_copies,
    unsigned const* lid_of_copies,
    unsigned const* msg_of_uses,
    unsigned const* rank_of_msgs)
{
  unsigned* lid_of_uses = LOOP_MALLOC(unsigned, nuses);
  for (unsigned i = 0; i < nowners; ++i) {
    unsigned fu = uses_of_owners_offsets[i];
    unsigned eu = uses_of_owners_offsets[i + 1];
    unsigned fc = copies_of_owners_offsets[i];
    unsigned ec = copies_of_owners_offsets[i + 1];
    for (unsigned j = fu; j < eu; ++j) {
      unsigned msg = msg_of_uses[j];
      unsigned rank = rank_of_msgs[msg];
      unsigned lid = INVALID;
      for (unsigned k = fc; k < ec; ++k)
        if (rank == rank_of_copies[k]) {
          lid = lid_of_copies[k];
          break;
        }
      assert(lid != INVALID);
      lid_of_uses[j] = lid;
    }
  }
  return lid_of_uses;
}

static struct mesh* migrate_element_topology(
    struct mesh* m,
    unsigned nelems_recvd,
    unsigned const* recvd_elem_ranks,
    unsigned const* recvd_elem_ids,
    /* this exchanger sends from new elements to old elements */
    struct exchanger** p_elem_pull,
    /* this exchanger sends from old vertex owners to new vertex copies */
    struct exchanger** p_vcopy_push)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  unsigned nverts_per_elem = the_down_degrees[dim][0];
  unsigned* or_of_vus_of_es;
  unsigned* oi_of_vus_of_es;
  get_vert_use_owners_of_elems(m, &or_of_vus_of_es, &oi_of_vus_of_es);
  struct exchanger* elem_pull = new_exchanger(nelems_recvd, recvd_elem_ranks);
  set_exchanger_dests(elem_pull, nelems, recvd_elem_ids);
  /* own rank of vertex uses of received elements */
  unsigned* or_of_vus_of_res = exchange_uints(elem_pull, nverts_per_elem,
      or_of_vus_of_es, EX_REV, EX_ITEM);
  loop_free(or_of_vus_of_es);
  unsigned* oi_of_vus_of_res = exchange_uints(elem_pull, nverts_per_elem,
      oi_of_vus_of_es, EX_REV, EX_ITEM);
  loop_free(oi_of_vus_of_es);
  unsigned nverts = mesh_count(m, 0);
  /* received vertex uses to old owner vertices exchanger */
  struct exchanger* us_to_own_ex = new_exchanger(nelems_recvd * nverts_per_elem,
      or_of_vus_of_res);
  loop_free(or_of_vus_of_res);
  set_exchanger_dests(us_to_own_ex, nverts, oi_of_vus_of_res);
  loop_free(oi_of_vus_of_res);
  unsigned* copies_of_owners_offsets;
  unsigned* rank_of_copies;
  unsigned const* rus_of_own_offsets =
    us_to_own_ex->items_of_roots_offsets[EX_REV];
  get_unique_ranks_of_owners(nverts, rus_of_own_offsets,
      us_to_own_ex->msg_of_items[EX_REV], us_to_own_ex->ranks[EX_REV],
      &copies_of_owners_offsets, &rank_of_copies);
  unsigned ncopies = copies_of_owners_offsets[nverts];
  struct exchanger* copies_push = new_exchanger(ncopies, rank_of_copies);
  set_exchanger_srcs(copies_push, nverts, copies_of_owners_offsets);
  unsigned nverts_recvd = copies_push->nitems[EX_REV];
  unsigned* lids = uints_linear(nverts_recvd);
  unsigned* lid_of_copies = exchange_uints(copies_push, 1, lids,
      EX_REV, EX_ITEM);
  loop_free(lids);
  unsigned nuses_recvd = us_to_own_ex->nitems[EX_REV];
  unsigned* lid_of_rvus = use_lids_from_copy_lids(
      nuses_recvd, nverts, rus_of_own_offsets,
      copies_of_owners_offsets, rank_of_copies,
      lid_of_copies, us_to_own_ex->msg_of_items[EX_REV],
      us_to_own_ex->ranks[EX_REV]);
  loop_free(copies_of_owners_offsets);
  loop_free(rank_of_copies);
  loop_free(lid_of_copies);
  unsigned* verts_of_elems = exchange_uints(us_to_own_ex, 1, lid_of_rvus,
      EX_REV, EX_ITEM);
  loop_free(lid_of_rvus);
  free_exchanger(us_to_own_ex);
  struct mesh* m_out = new_mesh(dim);
  mesh_set_ents(m_out, 0, nverts_recvd, 0);
  mesh_set_ents(m_out, dim, nelems_recvd, verts_of_elems);
  *p_elem_pull = elem_pull;
  *p_vcopy_push = copies_push;
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
      EX_REV, EX_ITEM);
  free_exchanger(elem_pull);
  free_mesh(*p_m);
  *p_m = m_out;
}
