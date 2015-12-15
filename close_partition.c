#include "close_partition.h"

#include "arrays.h"
#include "exchanger.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"

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

/* for each entity, the owner copy in the old mesh
   has just received a list of ranks that will require
   a copy in the new mesh, with duplicate ranks.
   this function and the above helper simply remove
   the duplicates and output, for each owner copy in
   the old mesh, the list of unique ranks that
   will have a copy of that entity in the new mesh */

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

/* we here have two entity dimensions A and B.
   the input is: for each A entity copy in this rank
   in the new mesh, the owners (in the old mesh)
   of the B entities it is adjacent to.
   This algorithm brings those adjacent entities
   into this rank, resulting in output arrays
   describing the B entity copies that will be
   in this rank, identified by their owners in
   the old mesh.
   the main example uses:
     given where elements go, determine where vertices go
     given where vertices go, determine where elements go (ghosting) */

void close_partition_exchangers(
    unsigned nacopies,
    unsigned nbowners,
    unsigned const* buses_by_acopies_offsets,
    unsigned const* buse_own_ranks,
    unsigned const* buse_own_ids,
    struct exchanger** p_buse_to_own,
    struct exchanger** p_bown_to_copy)
{
  unsigned nbuses = buses_by_acopies_offsets[nacopies];
  struct exchanger* buse_to_own = new_exchanger(nbuses, buse_own_ranks);
  set_exchanger_dests(buse_to_own, nbowners, buse_own_ids);
  unsigned* copies_of_owners_offsets;
  unsigned* rank_of_copies;
  get_unique_ranks_of_owners(nbowners,
      buse_to_own->items_of_roots_offsets[EX_REV],
      buse_to_own->msg_of_items[EX_REV],
      buse_to_own->ranks[EX_REV],
      &copies_of_owners_offsets,
      &rank_of_copies);
  unsigned nbcopies_own = copies_of_owners_offsets[nbowners];
  struct exchanger* bown_to_copy = new_exchanger(nbcopies_own,
      rank_of_copies);
  loop_free(rank_of_copies);
  set_exchanger_srcs(bown_to_copy, nbowners, copies_of_owners_offsets);
  loop_free(copies_of_owners_offsets);
  *p_buse_to_own = buse_to_own;
  *p_bown_to_copy = bown_to_copy;
}

void close_partition(
    unsigned nacopies,
    unsigned nbowners,
    unsigned const* buses_by_acopies_offsets,
    unsigned const* buse_own_ranks,
    unsigned const* buse_own_ids,
    unsigned* p_nbcopies,
    unsigned** p_bcopy_own_ranks,
    unsigned** p_bcopy_own_ids)
{
  struct exchanger* buse_to_own;
  struct exchanger* bown_to_copy;
  close_partition_exchangers(nacopies, nbowners, buses_by_acopies_offsets,
      buse_own_ranks, buse_own_ids, &buse_to_own, &bown_to_copy);
  free_exchanger(buse_to_own);
  unsigned nbcopies = bown_to_copy->nitems[EX_REV];
  unsigned* bcopy_own_ranks = LOOP_MALLOC(unsigned, nbcopies);
  for (unsigned i = 0; i < nbcopies; ++i)
    bcopy_own_ranks[i] = bown_to_copy->ranks[EX_REV][
      bown_to_copy->msg_of_items[EX_REV][i]];
  unsigned* lids = uints_linear(nbowners, 1);
  unsigned* bcopy_own_ids = exchange_uints(bown_to_copy, 1, lids,
      EX_FOR, EX_ROOT);
  loop_free(lids);
  free_exchanger(bown_to_copy);
  *p_nbcopies = nbcopies;
  *p_bcopy_own_ranks = bcopy_own_ranks;
  *p_bcopy_own_ids = bcopy_own_ids;
}

/* for each high enitty in the mesh, for each low entity it uses,
   return the owner of that used entity.
   the result is of size (nhighs * nlows_per_high) */

void get_down_use_owners(
    struct mesh* m,
    unsigned high_dim,
    unsigned low_dim,
    unsigned** p_use_own_ranks,
    unsigned** p_use_own_ids,
    unsigned** p_uses_of_highs_offsets)
{
  unsigned nhighs = mesh_count(m, high_dim);
  unsigned nlows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned const* lows_of_highs = mesh_ask_down(m, high_dim, low_dim);
  unsigned const* low_own_ranks = mesh_ask_own_ranks(m, low_dim);
  unsigned const* low_own_ids = mesh_ask_own_ids(m, low_dim);
  unsigned* use_own_ranks = LOOP_MALLOC(unsigned, nhighs * nlows_per_high);
  unsigned* use_own_ids = LOOP_MALLOC(unsigned, nhighs * nlows_per_high);
  for (unsigned i = 0; i < nhighs; ++i) {
    for (unsigned j = 0; j < nlows_per_high; ++j) {
      unsigned low = lows_of_highs[i * nlows_per_high + j];
      use_own_ranks[i * nlows_per_high + j] = low_own_ranks[low];
      use_own_ids[i * nlows_per_high + j] = low_own_ids[low];
    }
  }
  *p_use_own_ranks = use_own_ranks;
  *p_use_own_ids = use_own_ids;
  *p_uses_of_highs_offsets = uints_linear(nhighs, nlows_per_high);
}

/* this is supposed to be a simple operation.
   for each entity copy in the old mesh, we have a set of entities
   it uses, identified by their owners, and we
   simply want to communicate that set to all
   copies of that entity in the new mesh.
   unfortunately, the whole exchanger system is built for equally-sized
   entries (via the "width" argument), and in the cases of
   vertices "using" elements we have variable-size data.
   so this function is mostly a custom uints_expand for that
   possibly-variable-sized data.
   this also involves making another exchanger */

void pull_use_owners(
    struct exchanger* pull,
    unsigned const* use_own_ranks_in,
    unsigned const* use_own_ids_in,
    unsigned const* offsets_in,
    unsigned** p_use_own_ranks_out,
    unsigned** p_use_own_ids_out,
    unsigned** p_offsets_out)
{
  unsigned nents_out = pull->nitems[EX_FOR];
  unsigned nents_in = pull->nroots[EX_REV];
  unsigned const* out_of_in_offsets = pull->items_of_roots_offsets[EX_REV];
  unsigned* nout_of_in = uints_unscan(out_of_in_offsets, nents_in);
  unsigned* nuses_of_in = uints_unscan(offsets_in, nents_in);
  unsigned* prod_counts = LOOP_MALLOC(unsigned, nents_in);
  for (unsigned i = 0; i < nents_in; ++i)
    prod_counts[i] = nout_of_in[i] * nuses_of_in[i];
  loop_free(nout_of_in);
  loop_free(nuses_of_in);
  unsigned* prod_offsets = uints_exscan(prod_counts, nents_in);
  loop_free(prod_counts);
  unsigned* lids_out = uints_linear(nents_out, 1);
  unsigned* lids_in = exchange_uints(pull, 1, lids_out, EX_FOR, EX_ITEM);
  loop_free(lids_out);
  unsigned nprod = prod_offsets[nents_in];
  unsigned* prod_dest_ranks = LOOP_MALLOC(unsigned, nprod);
  unsigned* prod_dest_ids = LOOP_MALLOC(unsigned, nprod);
  unsigned* prod_use_ranks = LOOP_MALLOC(unsigned, nprod);
  unsigned* prod_use_ids = LOOP_MALLOC(unsigned, nprod);
  for (unsigned i = 0; i < nents_in; ++i) {
    unsigned fu = offsets_in[i];
    unsigned eu = offsets_in[i + 1];
    unsigned fo = out_of_in_offsets[i];
    unsigned eo = out_of_in_offsets[i + 1];
    unsigned l = prod_offsets[i];
    for (unsigned j = fu; j < eu; ++j)
      for (unsigned k = fo; k < eo; ++k) {
        prod_dest_ranks[l] = pull->ranks[EX_REV][
          pull->msg_of_items[EX_REV][k]];
        prod_dest_ids[l] = lids_in[k];
        prod_use_ranks[l] = use_own_ranks_in[j];
        prod_use_ids[l] = use_own_ids_in[j];
        ++l;
      }
  }
  loop_free(lids_in);
  loop_free(prod_offsets);
  struct exchanger* prod_push = new_exchanger(nprod, prod_dest_ranks);
  loop_free(prod_dest_ranks);
  set_exchanger_dests(prod_push, nents_out, prod_dest_ids);
  loop_free(prod_dest_ids);
  *p_use_own_ranks_out = exchange_uints(prod_push, 1, prod_use_ranks,
      EX_FOR, EX_ITEM);
  loop_free(prod_use_ranks);
  *p_use_own_ids_out = exchange_uints(prod_push, 1, prod_use_ids,
      EX_FOR, EX_ITEM);
  loop_free(prod_use_ids);
  *p_offsets_out = uints_copy(prod_push->items_of_roots_offsets[EX_REV],
      nents_out + 1);
  free_exchanger(prod_push);
}
