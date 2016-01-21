#include "parallel_modify.h"

#include <assert.h>

#include "arrays.h"
#include "comm.h"
#include "exchanger.h"
#include "ints.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"

/* given a ghosted mesh with an independent set
   of some key entities tagged as "indset",
   this function gives ownership of the elements
   adjacent to a key to the rank that owns the key */

void set_own_ranks_by_indset(
    struct mesh* m,
    unsigned key_dim)
{
  if (!mesh_is_parallel(m))
    return;
  assert(mesh_ghost_layers(m) == 1);
  unsigned const* indset = mesh_find_tag(m, key_dim, "indset")->d.u32;
  unsigned const* key_owners = mesh_ask_own_ranks(m, key_dim);
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned* elem_owners = uints_copy(
      mesh_ask_own_ranks(m, elem_dim), nelems);
  unsigned const* keys_of_elems = mesh_ask_down(m, elem_dim, key_dim);
  unsigned keys_per_elem = the_down_degrees[elem_dim][key_dim];
  for (unsigned i = 0; i < nelems; ++i) {
    for (unsigned j = 0; j < keys_per_elem; ++j) {
      unsigned key = keys_of_elems[i * keys_per_elem + j];
      if (indset[key])
        elem_owners[i] = key_owners[key];
    }
  }
  mesh_set_own_ranks(m, elem_dim, elem_owners);
}

void inherit_globals(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned ent_dim,
    unsigned const* offset_of_same_ents)
{
  /* this call ensures that the global numbers in the old mesh
     are organized such that the entities owned by the
     same rank are numbered consecutively */
  mesh_global_renumber(m_in, ent_dim);
  unsigned nents_in = mesh_count(m_in, ent_dim);
  unsigned nsame_ents = offset_of_same_ents[nents_in];
  unsigned* owned_ents_in = mesh_get_owned(m_in, ent_dim);
  unsigned nowned_ents_in = uints_sum(owned_ents_in, nents_in);
  unsigned* owned_same_ents = uints_expand(nents_in, 1,
      owned_ents_in, offset_of_same_ents);
  unsigned nowned_same_ents = uints_sum(
      owned_same_ents, nsame_ents);
  unsigned nents_out = mesh_count(m_out, ent_dim);
  unsigned nnew_ents = nents_out - nsame_ents;
  unsigned nowned_ents_out = nowned_same_ents + nnew_ents;
  unsigned long offset_out = comm_exscan_ulong(nowned_ents_out);
  unsigned long offset_in = comm_exscan_ulong(nowned_ents_in);
  long shift = ((long)offset_out) - ((long)offset_in);
  struct exchanger* ex = mesh_ask_exchanger(m_in, ent_dim);
  unsigned nowner_neighbors = ex->nmsgs[EX_REV];
  long* owner_neighbor_shifts = LOOP_MALLOC(long, nowner_neighbors);
  comm_sync_long(ex->comms[EX_FOR], shift, owner_neighbor_shifts);
  unsigned long* shifted_globals_in = ulongs_copy(
      mesh_ask_global(m_in, ent_dim), nents_in);
  for (unsigned i = 0; i < nents_in; ++i)
    shifted_globals_in[i] = (unsigned long)(
        ((long)(shifted_globals_in[i])) +
        owner_neighbor_shifts[ex->msg_of_items[EX_REV][i]]);
  unsigned long* globals_out = LOOP_MALLOC(unsigned long, nents_out);
  ulongs_expand_into(nents_in, 1, shifted_globals_in, offset_of_same_ents,
      globals_out);
  for (unsigned i = nsame_ents; i < nents_out; ++i)
    globals_out[i] = offset_out + i;
  mesh_set_global(m_out, ent_dim, globals_out);
}
