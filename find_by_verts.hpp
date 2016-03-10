#ifndef FIND_BY_VERTS_HPP
#define FIND_BY_VERTS_HPP

#include "loop.hpp"

/* strong optimization here due to the
   frequent use of these operations.
   we assume we are searching for edges
   or triangles, and have special code for each.
*/

LOOP_INOUT static inline unsigned find_tri_up(
    unsigned const* verts_wanted,
    unsigned const* verts_of_tris,
    unsigned const* tris_of_verts_offsets,
    unsigned const* tris_of_verts,
    unsigned const* tris_of_verts_directions)
{
  unsigned v0 = verts_wanted[0];
  unsigned f = tris_of_verts_offsets[v0];
  unsigned e = tris_of_verts_offsets[v0 + 1];
  for (unsigned i = f; i < e; ++i) {
    unsigned tri = tris_of_verts[i];
    unsigned dir = tris_of_verts_directions[i];
    unsigned const* verts_of_tri = verts_of_tris + tri * 3;
    if (((verts_of_tri[(dir + 1) % 3] == verts_wanted[1]) &&
         (verts_of_tri[(dir + 2) % 3] == verts_wanted[2])) ||
        ((verts_of_tri[(dir + 1) % 3] == verts_wanted[2]) &&
         (verts_of_tri[(dir + 2) % 3] == verts_wanted[1])))
      return tri;
  }
  LOOP_NORETURN(0);
}

LOOP_INOUT static inline unsigned find_edge_up(
    unsigned const* verts_wanted,
    unsigned const* verts_of_edges,
    unsigned const* edges_of_verts_offsets,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_directions)
{
  unsigned v0 = verts_wanted[0];
  unsigned f = edges_of_verts_offsets[v0];
  unsigned e = edges_of_verts_offsets[v0 + 1];
  for (unsigned i = f; i < e; ++i) {
    unsigned edge = edges_of_verts[i];
    unsigned dir = edges_of_verts_directions[i];
    unsigned const* verts_of_edge = verts_of_edges + edge * 2;
    if (verts_of_edge[1 - dir] == verts_wanted[1])
      return edge;
  }
  LOOP_NORETURN(0);
}

LOOP_INOUT static inline unsigned find_by_verts(
    unsigned ent_dim,
    unsigned const* verts_wanted,
    unsigned const* verts_of_ents,
    unsigned const* ents_of_verts_offsets,
    unsigned const* ents_of_verts,
    unsigned const* ents_of_verts_directions)
{
  switch (ent_dim) {
    case 0:
      return verts_wanted[0];
    case 1:
      return find_edge_up(verts_wanted, verts_of_ents,
          ents_of_verts_offsets, ents_of_verts, ents_of_verts_directions);
    case 2:
      return find_tri_up(verts_wanted, verts_of_ents,
          ents_of_verts_offsets, ents_of_verts, ents_of_verts_directions);
  };
  LOOP_NORETURN(0);
}

#endif
