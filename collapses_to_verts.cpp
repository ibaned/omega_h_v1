#include "collapses_to_verts.hpp"

#include <cassert>

#include "collapse_codes.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "tables.hpp"

namespace omega_h {

LOOP_KERNEL(collapse_to_vert,
    unsigned const* edges_of_verts_offsets,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_directions,
    unsigned const* col_codes,
    double const* col_quals_of_edges,
    unsigned* candidates,
    double* col_qual_of_verts)

  unsigned first_use = edges_of_verts_offsets[i];
  unsigned end_use = edges_of_verts_offsets[i + 1];
  double maxq = 0;
  unsigned is_collapsing = 0;
  for (unsigned j = first_use; j < end_use; ++j) {
    unsigned edge = edges_of_verts[j];
    unsigned direction = edges_of_verts_directions[j];
    if (!collapses(col_codes[edge], direction))
      continue;
    double q = col_quals_of_edges[edge * 2 + direction];
    if (!is_collapsing) {
      is_collapsing = 1;
      maxq = q;
    } else if (q > maxq) {
      maxq = q;
    }
  }
  candidates[i] = is_collapsing;
  if (is_collapsing)
    col_qual_of_verts[i] = maxq;
}

void valid_collapses_to_verts(struct mesh* m)
{
  unsigned const* col_codes = mesh_find_tag(m, 1, "col_codes")->d.u32;
  double const* col_quals_of_edges =
    mesh_find_tag(m, 1, "col_quals")->d.f64;
  unsigned nverts = mesh_count(m, 0);
  unsigned const* edges_of_verts_offsets =
    mesh_ask_up(m, 0, 1)->offsets;
  unsigned const* edges_of_verts =
    mesh_ask_up(m, 0, 1)->adj;
  unsigned const* edges_of_verts_directions =
    mesh_ask_up(m, 0, 1)->directions;
  unsigned* candidates = LOOP_MALLOC(unsigned, nverts);
  double* col_qual_of_verts = LOOP_MALLOC(double, nverts);
  LOOP_EXEC(collapse_to_vert, nverts,
      edges_of_verts_offsets,
      edges_of_verts,
      edges_of_verts_directions,
      col_codes,
      col_quals_of_edges,
      candidates,
      col_qual_of_verts);
  mesh_conform_array(m, 0, 1, &candidates);
  mesh_conform_array(m, 0, 1, &col_qual_of_verts);
  mesh_add_tag(m, 0, TAG_U32, "candidates", 1, candidates);
  mesh_add_tag(m, 0, TAG_F64, "col_qual", 1, col_qual_of_verts);
}

LOOP_KERNEL(vertex_dest,
    unsigned const* edges_of_verts_offsets,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_directions,
    unsigned const* verts_of_edges,
    unsigned const* indset,
    double const* col_qual_of_verts,
    unsigned const* col_codes,
    double const* col_quals_of_edges,
    unsigned* gen_vert_of_verts)
  if (!indset[i])
    return;
  unsigned first_use = edges_of_verts_offsets[i];
  unsigned end_use = edges_of_verts_offsets[i + 1];
  double maxq = col_qual_of_verts[i];
  if (maxq == 1e10)
    maxq = -1e10;
  unsigned gen_vert = INVALID;
  for (unsigned j = first_use; j < end_use; ++j) {
    unsigned edge = edges_of_verts[j];
    unsigned direction = edges_of_verts_directions[j];
    if (!collapses(col_codes[edge], direction))
      continue;
    double q = col_quals_of_edges[edge * 2 + direction];
    if (q == 1e10)
      q = -1e10;
    if (q == maxq) {
      gen_vert = verts_of_edges[edge * 2 + (1 - direction)];
      break;
    }
  }
  assert(gen_vert != INVALID);
  gen_vert_of_verts[i] = gen_vert;
}

unsigned* collapsing_vertex_destinations(struct mesh* m)
{
  unsigned const* col_codes = mesh_find_tag(m, 1, "col_codes")->d.u32;
  double const* col_quals_of_edges =
    mesh_find_tag(m, 1, "col_quals")->d.f64;
  unsigned const* indset = mesh_find_tag(m, 0, "indset")->d.u32;
  double const* col_qual_of_verts =
    mesh_find_tag(m, 0, "col_qual")->d.f64;
  unsigned nverts = mesh_count(m, 0);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* edges_of_verts_offsets =
    mesh_ask_up(m, 0, 1)->offsets;
  unsigned const* edges_of_verts =
    mesh_ask_up(m, 0, 1)->adj;
  unsigned const* edges_of_verts_directions =
    mesh_ask_up(m, 0, 1)->directions;
  unsigned* gen_vert_of_verts = LOOP_MALLOC(unsigned, nverts);
  LOOP_EXEC(vertex_dest, nverts,
      edges_of_verts_offsets,
      edges_of_verts,
      edges_of_verts_directions,
      verts_of_edges,
      indset,
      col_qual_of_verts,
      col_codes,
      col_quals_of_edges,
      gen_vert_of_verts);
  return gen_vert_of_verts;
}

}
