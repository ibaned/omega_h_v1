#include "concat.h"

#include <string.h>

#include "loop.h"
#include "subset.h"
#include "tables.h"

static void* generic_concat(
    unsigned typesize,
    unsigned width,
    void const* a,
    unsigned na,
    void const* b,
    unsigned nb)
{
  unsigned a_bytes = na * width * typesize;
  unsigned b_bytes = nb * width * typesize;
  unsigned total_bytes = a_bytes + b_bytes;
  char* out = LOOP_MALLOC(char, total_bytes);
  memcpy(out, a, a_bytes);
  memcpy(out + a_bytes, b, b_bytes);
  return out;
}

unsigned* concat_uints(
    unsigned width,
    unsigned const* a,
    unsigned na,
    unsigned const* b,
    unsigned nb)
{
  return (unsigned*) generic_concat(sizeof(unsigned), width, a, na, b, nb);
}

double* concat_doubles(
    unsigned width,
    double const* a,
    unsigned na,
    double const* b,
    unsigned nb)
{
  return (double*) generic_concat(sizeof(double), width, a, na, b, nb);
}

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
