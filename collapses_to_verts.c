#include "collapses_to_verts.h"

#include "collapse_codes.h"
#include "loop.h"
#include "tables.h"

void valid_collapses_to_verts(
    struct mesh* m,
    double const* col_quals_of_edges,
    unsigned** candidates_out,
    unsigned** gen_vert_of_verts_out,
    double** col_qual_of_verts_out)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* edges_of_verts_offsets =
    mesh_ask_up(m, 0, 1)->offsets;
  unsigned const* edges_of_verts =
    mesh_ask_up(m, 0, 1)->adj;
  unsigned const* edges_of_verts_directions =
    mesh_ask_up(m, 0, 1)->directions;
  unsigned const* col_codes = mesh_find_tag(m, 1, "col_codes");
  unsigned* candidates = LOOP_MALLOC(unsigned, nverts);
  unsigned* gen_vert_of_verts = LOOP_MALLOC(unsigned, nverts);
  double* col_qual_of_verts = LOOP_MALLOC(double, nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_use = edges_of_verts_offsets[i];
    unsigned end_use = edges_of_verts_offsets[i + 1];
    double maxq = 0;
    unsigned gen_vert = INVALID;
    unsigned is_collapsing = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned edge = edges_of_verts[j];
      unsigned direction = edges_of_verts_directions[j];
      if (!collapses(col_codes[edge], direction))
        continue;
      double q = col_quals_of_edges[edge * 2 + direction];
      unsigned other_vert = verts_of_edges[edge * 2 + (1 - direction)];
      if (!is_collapsing) {
        is_collapsing = 1;
        maxq = q;
        gen_vert = other_vert;
      } else if (q > maxq) {
        maxq = q;
        gen_vert = other_vert;
      }
    }
    candidates[i] = is_collapsing;
    if (is_collapsing) {
      gen_vert_of_verts[i] = gen_vert;
      col_qual_of_verts[i] = maxq;
    }
  }
  mesh_free_tag(m, 1, "col_codes");
  *candidates_out = candidates;
  *gen_vert_of_verts_out = gen_vert_of_verts;
  *col_qual_of_verts_out = col_qual_of_verts;
}
