#ifndef FIND_BY_VERTS_H
#define FIND_BY_VERTS_H

#include "tables.h"

static inline unsigned find_by_verts(
    unsigned verts_per_ent,
    unsigned const* verts_of_target,
    unsigned const* verts_of_ents,
    unsigned const* ents_of_verts,
    unsigned const* ents_of_verts_offsets)
{
  unsigned vert = verts_of_target[0];
  for (unsigned j = ents_of_verts_offsets[vert];
      j < ents_of_verts_offsets[vert + 1]; ++j) {
    unsigned ent = ents_of_verts[j];
    unsigned const* verts_of_ent = verts_of_ents + ent * verts_per_ent;
    unsigned nmatches = 0;
    for (unsigned k = 0; k < verts_per_ent; ++k)
      for (unsigned l = 0; l < verts_per_ent; ++l)
        if (verts_of_ent[k] == verts_of_target[l])
          ++nmatches;
    /* found a match, has the same vertices */
    if (nmatches == verts_per_ent)
      return ent;
  }
  return INVALID;
}

#endif
