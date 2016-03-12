#include "star.hpp"

#include <cassert>

#include "arrays.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tables.hpp"

static LOOP_IN unsigned get_ent_star_general(
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned lows_per_high,
    unsigned low,
    unsigned* star)
{
  unsigned first_up = highs_of_lows_offsets[low];
  unsigned last_up = highs_of_lows_offsets[low + 1];
  unsigned size = 0;
  for (unsigned i = first_up; i < last_up; ++i) {
    unsigned high = highs_of_lows[i];
    for (unsigned j = 0; j < lows_per_high; ++j) {
      unsigned star_low = lows_of_highs[high * lows_per_high + j];
      if (star_low == low)
        continue;
      size = add_unique(star, size, star_low);
    }
  }
  return size;
}

LOOP_KERNEL(count_general,
    unsigned *degrees,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned lows_per_high)

  unsigned star_buf[MAX_UP * MAX_DOWN];
  degrees[i] = get_ent_star_general(
      highs_of_lows_offsets,
      highs_of_lows,
      lows_of_highs,
      lows_per_high,
      i,
      star_buf);
}

LOOP_KERNEL(fill_general,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned lows_per_high,
    unsigned* star,
    unsigned const* star_offsets)

  get_ent_star_general(
      highs_of_lows_offsets,
      highs_of_lows,
      lows_of_highs,
      lows_per_high,
      i,
      star + star_offsets[i]);
}

static void get_star_general(
    unsigned low_dim,
    unsigned high_dim,
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned** star_offsets_out,
    unsigned** star_out)
{
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned* degrees = filled_array<unsigned>(nlows, 0);
  LOOP_EXEC(count_general,
    nlows,
    degrees,
    highs_of_lows_offsets,
    highs_of_lows,
    lows_of_highs,
    lows_per_high);
  unsigned* star_offsets = uints_exscan(degrees, nlows);
  loop_free(degrees);
  unsigned sum_degrees = array_at(star_offsets, nlows);
  unsigned* star = LOOP_MALLOC(unsigned, sum_degrees);
  LOOP_EXEC(fill_general,
    nlows,
    highs_of_lows_offsets,
    highs_of_lows,
    lows_of_highs,
    lows_per_high,
    star,
    star_offsets);
  *star_offsets_out = star_offsets;
  *star_out = star;
}

static void mesh_get_star_general(
    struct mesh* m,
    unsigned low_dim,
    unsigned high_dim,
    unsigned** p_star_offsets,
    unsigned** p_star)
{
  get_star_general(low_dim, high_dim, mesh_count(m, low_dim),
      mesh_ask_up(m, low_dim, high_dim)->offsets,
      mesh_ask_up(m, low_dim, high_dim)->adj,
      mesh_ask_down(m, high_dim, low_dim),
      p_star_offsets, p_star);
}

LOOP_KERNEL(vertex_edge_kern,
    unsigned const* edges_of_verts_offsets,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_directions,
    unsigned const* verts_of_edges,
    unsigned* star)

  unsigned f = edges_of_verts_offsets[i];
  unsigned e = edges_of_verts_offsets[i + 1];
  for (unsigned j = f; j < e; ++j) {
    unsigned edge = edges_of_verts[j];
    unsigned dir = edges_of_verts_directions[j];
    star[j] = verts_of_edges[edge * 2 + (1 - dir)];
  }
}

static void get_vertex_edge_star(
    struct mesh* m,
    unsigned** p_star_offsets,
    unsigned** p_star)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned const* verts_of_edges =
    mesh_ask_down(m, 1, 0);
  unsigned const* edges_of_verts_offsets =
    mesh_ask_up(m, 0, 1)->offsets;
  unsigned const* edges_of_verts =
    mesh_ask_up(m, 0, 1)->adj;
  unsigned const* edges_of_verts_directions =
    mesh_ask_up(m, 0, 1)->directions;
  *p_star_offsets = copy_array(
      edges_of_verts_offsets, nverts + 1);
  unsigned nadj = array_at(edges_of_verts_offsets, nverts);
  unsigned* star = LOOP_MALLOC(unsigned, nadj);
  LOOP_EXEC(vertex_edge_kern, nverts,
      edges_of_verts_offsets,
      edges_of_verts,
      edges_of_verts_directions,
      verts_of_edges,
      star);
  *p_star = star;
}

LOOP_KERNEL(edge_tri_kern_1,
    unsigned const* tris_of_edges_offsets,
    unsigned* star_offsets)

  star_offsets[i] = tris_of_edges_offsets[i] * 2;
}

LOOP_KERNEL(edge_tri_kern_2,
    unsigned const* tris_of_edges_offsets,
    unsigned const* tris_of_edges,
    unsigned const* tris_of_edges_directions,
    unsigned const* edges_of_tris,
    unsigned* star)

  unsigned f = tris_of_edges_offsets[i];
  unsigned e = tris_of_edges_offsets[i + 1];
  for (unsigned j = f; j < e; ++j) {
    unsigned tri = tris_of_edges[j];
    unsigned dir = tris_of_edges_directions[j];
    star[j * 2 + 0] = edges_of_tris[tri * 3 + ((dir + 1) % 3)];
    star[j * 2 + 1] = edges_of_tris[tri * 3 + ((dir + 2) % 3)];
  }
}

static void get_edge_triangle_star(
    struct mesh* m,
    unsigned** p_star_offsets,
    unsigned** p_star)
{
  unsigned nedges = mesh_count(m, 1);
  unsigned const* edges_of_tris =
    mesh_ask_down(m, 2, 1);
  unsigned const* tris_of_edges_offsets =
    mesh_ask_up(m, 1, 2)->offsets;
  unsigned const* tris_of_edges =
    mesh_ask_up(m, 1, 2)->adj;
  unsigned const* tris_of_edges_directions =
    mesh_ask_up(m, 1, 2)->directions;
  unsigned* star_offsets = LOOP_MALLOC(unsigned, nedges + 1);
  LOOP_EXEC(edge_tri_kern_1, nedges + 1,
      tris_of_edges_offsets, star_offsets);
  unsigned nadj = array_at(tris_of_edges_offsets, nedges) * 2;
  unsigned* star = LOOP_MALLOC(unsigned, nadj);
  LOOP_EXEC(edge_tri_kern_2, nedges,
      tris_of_edges_offsets,
      tris_of_edges,
      tris_of_edges_directions,
      edges_of_tris,
      star);
  *p_star_offsets = star_offsets;
  *p_star = star;
}

LOOP_KERNEL(edge_tet_kern_1,
    unsigned const* edge_tri_star_degrees,
    unsigned const* edge_tet_degrees,
    unsigned* star_degrees)

  star_degrees[i] = edge_tri_star_degrees[i] + edge_tet_degrees[i];
}

LOOP_KERNEL(edge_tet_kern_2,
    unsigned const* star_offsets,
    unsigned const* edge_tri_star_offsets,
    unsigned const* edge_tri_star,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* edges_of_tets,
    unsigned* star)

  unsigned const* tet_edge_opp_edges = the_opposite_orders[3][1];
  unsigned o = star_offsets[i];
  {
    unsigned f = edge_tri_star_offsets[i];
    unsigned e = edge_tri_star_offsets[i + 1];
    for (unsigned j = f; j < e; ++j)
      star[o++] = edge_tri_star[j];
  }
  {
    unsigned f = tets_of_edges_offsets[i];
    unsigned e = tets_of_edges_offsets[i + 1];
    for (unsigned j = f; j < e; ++j) {
      unsigned tet = tets_of_edges[j];
      unsigned dir = tets_of_edges_directions[j];
      unsigned tet_edge = tet_edge_opp_edges[dir];
      assert(edges_of_tets[tet * 6 + dir] == i);
      star[o++] = edges_of_tets[tet * 6 + tet_edge];
    }
  }
  assert(o == star_offsets[i + 1]);
}

static void get_edge_tet_star(
    struct mesh* m,
    unsigned** p_star_offsets,
    unsigned** p_star)
{
  unsigned* edge_tri_star_offsets;
  unsigned* edge_tri_star;
  get_edge_triangle_star(m, &edge_tri_star_offsets, &edge_tri_star);
  unsigned nedges = mesh_count(m, 1);
  unsigned const* edges_of_tets =
    mesh_ask_down(m, 3, 1);
  unsigned const* tets_of_edges_offsets =
    mesh_ask_up(m, 1, 3)->offsets;
  unsigned const* tets_of_edges =
    mesh_ask_up(m, 1, 3)->adj;
  unsigned const* tets_of_edges_directions =
    mesh_ask_up(m, 1, 3)->directions;
  unsigned* edge_tri_star_degrees = uints_unscan(
      edge_tri_star_offsets, nedges);
  unsigned* edge_tet_degrees = uints_unscan(
      tets_of_edges_offsets, nedges);
  unsigned* star_degrees = LOOP_MALLOC(unsigned, nedges);
  LOOP_EXEC(edge_tet_kern_1, nedges,
      edge_tri_star_degrees, edge_tet_degrees, star_degrees);
  loop_free(edge_tri_star_degrees);
  loop_free(edge_tet_degrees);
  unsigned* star_offsets = uints_exscan(star_degrees, nedges);
  loop_free(star_degrees);
  unsigned nadj = array_at(star_offsets, nedges);
  unsigned* star = filled_array<unsigned>(nadj, INVALID);
  LOOP_EXEC(edge_tet_kern_2, nedges,
      star_offsets, edge_tri_star_offsets, edge_tri_star,
      tets_of_edges_offsets, tets_of_edges, tets_of_edges_directions,
      edges_of_tets, star);
  loop_free(edge_tri_star_offsets);
  loop_free(edge_tri_star);
  *p_star_offsets = star_offsets;
  *p_star = star;
}

void mesh_get_star(
    struct mesh* m,
    unsigned low_dim,
    unsigned high_dim,
    unsigned** p_star_offsets,
    unsigned** p_star)
{
  if (low_dim == 0 &&
      high_dim == mesh_dim(m) &&
      mesh_dim(m) != 1 &&
      mesh_has_dim(m, 1))
    mesh_get_star(m, low_dim, 1, p_star_offsets, p_star);
  else if (low_dim == 0 && high_dim == 1)
    get_vertex_edge_star(m, p_star_offsets, p_star);
  else if (low_dim == 1 && high_dim == 2)
    get_edge_triangle_star(m, p_star_offsets, p_star);
  else if (low_dim == 1 && high_dim == 3 && mesh_has_dim(m, 2))
    get_edge_tet_star(m, p_star_offsets, p_star);
  else
    mesh_get_star_general(m, low_dim, high_dim, p_star_offsets, p_star);
}
