#include "concat_topology.h"

#include <string.h>

#include "arrays.h"
#include "loop.h"
#include "subset.h"
#include "tables.h"

void concat_verts_of_ents(
    unsigned ent_dim,
    unsigned nents,
    unsigned ngen_ents,
    unsigned const* verts_of_ents,
    unsigned const* offset_of_same_ents,
    unsigned const* verts_of_gen_ents,
    unsigned* nents_out,
    unsigned** verts_of_ents_out)
{
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned nsame_ents = offset_of_same_ents[nents];
  unsigned* verts_of_same_ents = uints_subset(nents, verts_per_ent,
      verts_of_ents, offset_of_same_ents);
  unsigned* out = concat_uints(verts_per_ent,
      verts_of_same_ents, nsame_ents,
      verts_of_gen_ents, ngen_ents);
  loop_free(verts_of_same_ents);
  *nents_out = nsame_ents + ngen_ents;
  *verts_of_ents_out = out;
}
