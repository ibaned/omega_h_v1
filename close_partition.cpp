#include "close_partition.hpp"

#include "arrays.hpp"
#include "exchanger.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "tables.hpp"

namespace omega_h {

LOOP_IN static unsigned get_unique_ranks_of_owner(
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

LOOP_KERNEL(unique_ranks_1,
    unsigned const* uses_of_owners_offsets,
    unsigned const* msg_of_uses,
    unsigned const* rank_of_msgs,
    unsigned* ncopies_of_owners,
    unsigned* scratch)
  ncopies_of_owners[i] = get_unique_ranks_of_owner(
      uses_of_owners_offsets, msg_of_uses, rank_of_msgs,
      i, scratch + uses_of_owners_offsets[i]);
}

LOOP_KERNEL(unique_ranks_2,
    unsigned const* uses_of_owners_offsets,
    unsigned const* copies_of_owners_offsets,
    unsigned const* scratch,
    unsigned* rank_of_copies)
  unsigned fu = uses_of_owners_offsets[i];
  unsigned fc = copies_of_owners_offsets[i];
  unsigned ec = copies_of_owners_offsets[i + 1];
  unsigned nc = ec - fc;
  for (unsigned j = 0; j < nc; ++j)
    rank_of_copies[fc + j] = scratch[fu + j];
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
  unsigned nuses = array_at(uses_of_owners_offsets, nowners);
  unsigned* scratch = LOOP_MALLOC(unsigned, nuses);
  LOOP_EXEC(unique_ranks_1, nowners, uses_of_owners_offsets,
      msg_of_uses, rank_of_msgs, ncopies_of_owners, scratch);
  unsigned* copies_of_owners_offsets = uints_exscan(
      ncopies_of_owners, nowners);
  loop_free(ncopies_of_owners);
  unsigned ncopies = array_at(copies_of_owners_offsets, nowners);
  unsigned* rank_of_copies = LOOP_MALLOC(unsigned, ncopies);
  LOOP_EXEC(unique_ranks_2, nowners, uses_of_owners_offsets,
      copies_of_owners_offsets, scratch, rank_of_copies);
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

struct exchanger* close_partition_exchanger(
    struct exchanger* buse_to_own)
{
  unsigned nbowners = buse_to_own->nroots[EX_REV];
  unsigned* copies_of_owners_offsets;
  unsigned* rank_of_copies;
  get_unique_ranks_of_owners(nbowners,
      buse_to_own->items_of_roots_offsets[EX_REV],
      buse_to_own->msg_of_items[EX_REV],
      buse_to_own->ranks[EX_REV],
      &copies_of_owners_offsets,
      &rank_of_copies);
  unsigned nbcopies_own = array_at(copies_of_owners_offsets, nbowners);
  struct exchanger* bown_to_copy = new_exchanger(nbcopies_own,
      rank_of_copies);
  loop_free(rank_of_copies);
  set_exchanger_srcs(bown_to_copy, nbowners, copies_of_owners_offsets);
  loop_free(copies_of_owners_offsets);
  return bown_to_copy;
}

LOOP_KERNEL(msg_to_rank,
    unsigned const* msg_ranks,
    unsigned const* msg_of_items,
    unsigned* item_ranks)
  item_ranks[i] = msg_ranks[msg_of_items[i]];
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
  unsigned nuses = array_at(buses_by_acopies_offsets, nacopies);
  struct exchanger* buse_to_own = new_exchanger(nuses, buse_own_ranks);
  set_exchanger_dests(buse_to_own, nbowners, buse_own_ids);
  struct exchanger* bown_to_copy = close_partition_exchanger(buse_to_own);
  free_exchanger(buse_to_own);
  unsigned nbcopies = bown_to_copy->nitems[EX_REV];
  unsigned* bcopy_own_ranks = LOOP_MALLOC(unsigned, nbcopies);
  LOOP_EXEC(msg_to_rank, nbcopies,
      bown_to_copy->ranks[EX_REV],
      bown_to_copy->msg_of_items[EX_REV],
      bcopy_own_ranks);
  unsigned* lids = uints_linear(nbowners, 1);
  unsigned* bcopy_own_ids = exchange(bown_to_copy, 1, lids,
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

LOOP_KERNEL(down_owners,
    unsigned const* lows_of_highs,
    unsigned nlows_per_high,
    unsigned const* low_own_ranks,
    unsigned const* low_own_ids,
    unsigned* use_own_ranks,
    unsigned* use_own_ids)
  for (unsigned j = 0; j < nlows_per_high; ++j) {
    unsigned low = lows_of_highs[i * nlows_per_high + j];
    use_own_ranks[i * nlows_per_high + j] = low_own_ranks[low];
    use_own_ids[i * nlows_per_high + j] = low_own_ids[low];
  }
}

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
  LOOP_EXEC(down_owners, nhighs, lows_of_highs, nlows_per_high,
      low_own_ranks, low_own_ids,
      use_own_ranks, use_own_ids);
  *p_use_own_ranks = use_own_ranks;
  *p_use_own_ids = use_own_ids;
  *p_uses_of_highs_offsets = uints_linear(nhighs + 1, nlows_per_high);
}

LOOP_KERNEL(count_prods,
    unsigned const* nout_of_in,
    unsigned const* nuses_of_in,
    unsigned* prod_counts)
  prod_counts[i] = nout_of_in[i] * nuses_of_in[i];
}

LOOP_KERNEL(fill_prods,
    unsigned const* offsets_in,
    unsigned const* out_of_in_offsets,
    unsigned const* prod_offsets,
    unsigned const* ranks,
    unsigned const* msg_of_items,
    unsigned const* lids_in,
    unsigned const* use_own_ranks_in,
    unsigned const* use_own_ids_in,
    unsigned* prod_dest_ranks,
    unsigned* prod_dest_ids,
    unsigned* prod_use_ranks,
    unsigned* prod_use_ids)
  unsigned fu = offsets_in[i];
  unsigned eu = offsets_in[i + 1];
  unsigned fo = out_of_in_offsets[i];
  unsigned eo = out_of_in_offsets[i + 1];
  unsigned l = prod_offsets[i];
  for (unsigned j = fu; j < eu; ++j)
    for (unsigned k = fo; k < eo; ++k) {
      prod_dest_ranks[l] = ranks[msg_of_items[k]];
      prod_dest_ids[l] = lids_in[k];
      prod_use_ranks[l] = use_own_ranks_in[j];
      prod_use_ids[l] = use_own_ids_in[j];
      ++l;
    }
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

void push_use_owners(
    struct exchanger* push,
    unsigned const* use_own_ranks_in,
    unsigned const* use_own_ids_in,
    unsigned const* offsets_in,
    unsigned** p_use_own_ranks_out,
    unsigned** p_use_own_ids_out,
    unsigned** p_offsets_out)
{
  unsigned nents_out = push->nitems[EX_REV];
  unsigned nents_in = push->nroots[EX_FOR];
  unsigned const* out_of_in_offsets = push->items_of_roots_offsets[EX_FOR];
  unsigned* nout_of_in = uints_unscan(out_of_in_offsets, nents_in);
  unsigned* nuses_of_in = uints_unscan(offsets_in, nents_in);
  unsigned* prod_counts = LOOP_MALLOC(unsigned, nents_in);
  LOOP_EXEC(count_prods, nents_in, nout_of_in, nuses_of_in, prod_counts);
  loop_free(nout_of_in);
  loop_free(nuses_of_in);
  unsigned* prod_offsets = uints_exscan(prod_counts, nents_in);
  loop_free(prod_counts);
  unsigned* lids_out = uints_linear(nents_out, 1);
  unsigned* lids_in = exchange(push, 1, lids_out, EX_REV, EX_ITEM);
  loop_free(lids_out);
  unsigned nprod = array_at(prod_offsets, nents_in);
  unsigned* prod_dest_ranks = LOOP_MALLOC(unsigned, nprod);
  unsigned* prod_dest_ids = LOOP_MALLOC(unsigned, nprod);
  unsigned* prod_use_ranks = LOOP_MALLOC(unsigned, nprod);
  unsigned* prod_use_ids = LOOP_MALLOC(unsigned, nprod);
  LOOP_EXEC(fill_prods, nents_in,
      offsets_in,
      out_of_in_offsets,
      prod_offsets,
      /* FIXME: are these two arrays on device ?? */
      push->ranks[EX_FOR],
      push->msg_of_items[EX_FOR],
      lids_in,
      use_own_ranks_in,
      use_own_ids_in,
      prod_dest_ranks,
      prod_dest_ids,
      prod_use_ranks,
      prod_use_ids);
  loop_free(lids_in);
  loop_free(prod_offsets);
  struct exchanger* prod_push = new_exchanger(nprod, prod_dest_ranks);
  loop_free(prod_dest_ranks);
  set_exchanger_dests(prod_push, nents_out, prod_dest_ids);
  loop_free(prod_dest_ids);
  *p_use_own_ranks_out = exchange(prod_push, 1, prod_use_ranks,
      EX_FOR, EX_ITEM);
  loop_free(prod_use_ranks);
  *p_use_own_ids_out = exchange(prod_push, 1, prod_use_ids,
      EX_FOR, EX_ITEM);
  loop_free(prod_use_ids);
  *p_offsets_out = copy_array(prod_push->items_of_roots_offsets[EX_REV],
      nents_out + 1);
  free_exchanger(prod_push);
}

}
