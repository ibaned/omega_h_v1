#include "parallel_modify.h"

#include <assert.h>

#include "arrays.h"
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
