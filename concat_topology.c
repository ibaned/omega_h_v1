#include "concat_topology.h"

#include <string.h>

#include "arrays.h"
#include "loop.h"
#include "tables.h"

static void append(unsigned* out, unsigned* offset,
    unsigned const* in, unsigned n, unsigned width)
{
  unsigned nuints = n * width;
  unsigned nbytes = nuints * sizeof(unsigned);
  loop_memcpy(out + *offset, in, nbytes);
  *offset += nuints;
}

void concat_verts_of_ents(
    unsigned ent_dim,
    unsigned nents,
    unsigned const* verts_of_ents,
    unsigned const* offset_of_same_ents,
    unsigned const ngen_ents[4],
    unsigned* verts_of_gen_ents[4],
    unsigned* p_nents_out,
    unsigned** p_verts_of_ents_out)
{
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned nsame_ents = offset_of_same_ents[nents];
  unsigned* verts_of_same_ents = uints_expand(nents, verts_per_ent,
      verts_of_ents, offset_of_same_ents);
  unsigned nents_out = nsame_ents;
  for (unsigned i = 0; i < 4; ++i)
    nents_out += ngen_ents[i];
  unsigned* verts_of_ents_out = LOOP_MALLOC(unsigned,
      nents_out * verts_per_ent);
  unsigned offset = 0;
  append(verts_of_ents_out, &offset, verts_of_same_ents, nsame_ents,
      verts_per_ent);
  loop_free(verts_of_same_ents);
  for (unsigned i = 0; i < 4; ++i)
    append(verts_of_ents_out, &offset, verts_of_gen_ents[i], ngen_ents[i],
        verts_per_ent);
  *p_nents_out = nents_out;
  *p_verts_of_ents_out = verts_of_ents_out;
}

void concat_verts_of_elems(
    unsigned elem_dim,
    unsigned nelems,
    unsigned ngen_elems,
    unsigned const* verts_of_elems,
    unsigned const* offset_of_same_elems,
    unsigned* verts_of_gen_elems,
    unsigned* nelems_out,
    unsigned** verts_of_elems_out)
{
  unsigned ngen_ents[4] = {0,0,0,ngen_elems};
  unsigned* verts_of_gen_ents[4] =
  {0,0,0,verts_of_gen_elems};
  concat_verts_of_ents(elem_dim, nelems, verts_of_elems,
      offset_of_same_elems, ngen_ents, verts_of_gen_ents,
      nelems_out, verts_of_elems_out);
}
