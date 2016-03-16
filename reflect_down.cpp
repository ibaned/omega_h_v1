#include "reflect_down.hpp"

#include <cassert>

#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tables.hpp"

namespace omega_h {

typedef Kokkos::View<unsigned const*, Kokkos::MemoryUnmanaged> UCV;
typedef Kokkos::View<unsigned*, Kokkos::MemoryUnmanaged> UV;

LOOP_INOUT static inline unsigned find_tri_up(
    unsigned const* verts_wanted,
    UCV verts_of_tris,
    UCV tris_of_verts_offsets,
    UCV tris_of_verts,
    UCV tris_of_verts_directions)
{
  unsigned v0 = verts_wanted[0];
  unsigned f = tris_of_verts_offsets[v0];
  unsigned e = tris_of_verts_offsets[v0 + 1];
  for (unsigned i = f; i < e; ++i) {
    unsigned tri = tris_of_verts[i];
    unsigned dir = tris_of_verts_directions[i];
    unsigned const* verts_of_tri = &verts_of_tris[tri * 3];
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
    UCV verts_of_edges,
    UCV edges_of_verts_offsets,
    UCV edges_of_verts,
    UCV edges_of_verts_directions)
{
  unsigned v0 = verts_wanted[0];
  unsigned f = edges_of_verts_offsets[v0];
  unsigned e = edges_of_verts_offsets[v0 + 1];
  for (unsigned i = f; i < e; ++i) {
    unsigned edge = edges_of_verts[i];
    unsigned dir = edges_of_verts_directions[i];
    unsigned const* verts_of_edge = &verts_of_edges[edge * 2];
    if (verts_of_edge[1 - dir] == verts_wanted[1])
      return edge;
  }
  LOOP_NORETURN(0);
}

LOOP_INOUT static inline unsigned find_by_verts(
    unsigned ent_dim,
    unsigned const* verts_wanted,
    UCV verts_of_ents,
    UCV ents_of_verts_offsets,
    UCV ents_of_verts,
    UCV ents_of_verts_directions)
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

LOOP_KERNEL(reflect_down_entity,
    unsigned high_dim,
    unsigned low_dim,
    UCV verts_of_highs,
    UCV verts_of_lows,
    UCV lows_of_verts_offsets,
    UCV lows_of_verts,
    UCV lows_of_verts_directions,
    unsigned verts_per_high,
    unsigned lows_per_high,
    unsigned verts_per_low,
    UV lows_of_highs)

  unsigned const* const* high_verts_of_lows =
      the_canonical_orders[high_dim][low_dim][0];
  unsigned const* verts_of_high = &verts_of_highs[i * verts_per_high];
  unsigned* lows_of_high = &lows_of_highs[i * lows_per_high];
  for (unsigned j = 0; j < lows_per_high; ++j) {
    unsigned const* high_verts_of_low = high_verts_of_lows[j];
    unsigned verts_wanted[3];
    for (unsigned k = 0; k < verts_per_low; ++k)
      verts_wanted[k] = verts_of_high[high_verts_of_low[k]];
    lows_of_high[j] = find_by_verts(low_dim, verts_wanted, verts_of_lows,
        lows_of_verts_offsets, lows_of_verts, lows_of_verts_directions);
  }
}

static unsigned* reflect_down(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    UCV verts_of_highs,
    UCV verts_of_lows,
    UCV lows_of_verts_offsets,
    UCV lows_of_verts,
    UCV lows_of_verts_directions)
{
  unsigned verts_per_high = the_down_degrees[high_dim][0];
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned verts_per_low = the_down_degrees[low_dim][0];
  assert(verts_per_low == 2 || verts_per_low == 3);
  UV lows_of_highs(LOOP_MALLOC(unsigned, nhighs * lows_per_high),
                   nhighs * lows_per_high);
  LOOP_EXEC(reflect_down_entity, nhighs,
      high_dim,
      low_dim,
      verts_of_highs,
      verts_of_lows,
      lows_of_verts_offsets,
      lows_of_verts,
      lows_of_verts_directions,
      verts_per_high,
      lows_per_high,
      verts_per_low,
      lows_of_highs);
  return lows_of_highs.ptr_on_device();
}

unsigned* mesh_reflect_down(
    struct mesh* m,
    unsigned high_dim,
    unsigned low_dim)
{
  unsigned nhighs = mesh_count(m, high_dim);
  UCV verts_of_highs(mesh_ask_down(m, high_dim, 0),
      mesh_count(m, high_dim) * the_down_degrees[high_dim][0]);
  struct const_up* lows_of_verts = mesh_ask_up(m, 0, low_dim);
  UCV lows_of_verts_offsets(lows_of_verts->offsets,
      mesh_count(m, low_dim));
  UCV lows_of_verts_adj(lows_of_verts->adj,
      mesh_count(m, low_dim) * the_down_degrees[low_dim][0]);
  UCV lows_of_verts_directions(lows_of_verts->directions,
      mesh_count(m, low_dim) * the_down_degrees[low_dim][0]);
  UCV verts_of_lows(mesh_ask_down(m, low_dim, 0),
      mesh_count(m, low_dim) * the_down_degrees[low_dim][0]);
  return reflect_down(high_dim, low_dim, nhighs,
      verts_of_highs,
      verts_of_lows,
      lows_of_verts_offsets,
      lows_of_verts_adj,
      lows_of_verts_directions);
}

}
