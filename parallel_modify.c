#include "parallel_modify.h"

#include <assert.h>
#include <stdio.h>

#include "arrays.h"
#include "comm.h"
#include "exchanger.h"
#include "ints.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"
#include "vtk.h"

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
  unsigned nin = mesh_count(m_in, ent_dim);
  unsigned* owned = mesh_get_owned(m_in, ent_dim);
  unsigned* same = uints_unscan(offset_of_same_ents, nin);
  unsigned* owned_and_same = LOOP_MALLOC(unsigned, nin);
  for (unsigned i = 0; i < nin; ++i)
    owned_and_same[i] = owned[i] && same[i];
  loop_free(owned);
  loop_free(same);
  unsigned* owned_and_same_offsets = uints_exscan(owned_and_same, nin);
  unsigned nowned_and_same = owned_and_same_offsets[nin];
  unsigned nsame = offset_of_same_ents[nin];
  unsigned nout = mesh_count(m_out, ent_dim);
  unsigned nnew = nout - nsame;
  unsigned nowned_out = nowned_and_same + nnew;
  unsigned long offset_out = comm_exscan_ulong(nowned_out);
  unsigned long* new_globals_in = ulongs_filled(nin, ~((unsigned long) 0));
  for (unsigned i = 0; i < nin; ++i)
    if (owned_and_same[i])
      new_globals_in[i] = offset_out + owned_and_same_offsets[i];
  loop_free(owned_and_same);
  loop_free(owned_and_same_offsets);
  mesh_conform_ulongs(m_in, ent_dim, 1, &new_globals_in);
  unsigned long* new_globals_out = LOOP_MALLOC(unsigned long, nout);
  ulongs_expand_into(nin, 1, new_globals_in, offset_of_same_ents,
      new_globals_out);
  loop_free(new_globals_in);
  for (unsigned i = 0; i < nnew; ++i)
    new_globals_out[i + nsame] = offset_out + nowned_and_same + i;
  mesh_set_globals(m_out, ent_dim, new_globals_out);
}
