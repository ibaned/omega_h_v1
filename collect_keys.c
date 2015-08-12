#include "collect_keys.h"
#include "tables.h"
#include <stdlib.h>

unsigned* collect_keys(
    unsigned elem_dim,
    unsigned key_dim,
    unsigned nkeys,
    unsigned const* elems_of_keys_offsets,
    unsigned const* elems_of_keys,
    unsigned const* elems_of_keys_directions,
    unsigned const* bad_elems,
    unsigned const* key_of_elems)
{
  unsigned* candidates = malloc(sizeof(unsigned) * nkeys);
  /* handles the "double-split" behavior for sliver tets */
  unsigned opposites_too = (key_dim == get_opposite_dim(elem_dim, key_dim));
  unsigned const* elem_key_opp_key = the_opposite_orders[elem_dim][key_dim];
  for (unsigned i = 0; i < nkeys; ++i) {
    unsigned first_use = elems_of_keys_offsets[i];
    unsigned end_use = elems_of_keys_offsets[i + 1];
    unsigned is_candidate = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned elem = elems_of_keys[j];
      if (!bad_elems[elem])
        continue;
      unsigned direction = elems_of_keys_directions[j];
      if (direction == key_of_elems[elem])
        is_candidate = 1;
      if (!opposites_too)
        continue;
      unsigned opp = elem_key_opp_key[direction];
      if (opp == key_of_elems[elem])
        is_candidate = 1;
    }
    candidates[i] = is_candidate;
  }
  return candidates;
}
