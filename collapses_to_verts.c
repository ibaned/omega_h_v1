#include "collapses_to_verts.h"
#include "loop.h"  // for malloc
#include "tables.h"  // for INVALID

void collapses_to_verts(
    unsigned nverts,
    unsigned const* verts_of_edges,
    unsigned const* edges_of_verts_offsets,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_directions,
    unsigned const* col_codes,
    double const* col_quals_of_edges,
    unsigned** candidates_out,
    unsigned** gen_vert_of_verts_out,
    double** col_qual_of_verts_out)
{
  unsigned* candidates = loop_malloc(sizeof(unsigned) * nverts);
  unsigned* gen_vert_of_verts = loop_malloc(sizeof(unsigned) * nverts);
  double* col_qual_of_verts = loop_malloc(sizeof(double) * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_use = edges_of_verts_offsets[i];
    unsigned end_use = edges_of_verts_offsets[i + 1];
    double maxq = 0;
    unsigned gen_vert = INVALID;
    unsigned is_collapsing = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned edge = edges_of_verts[j];
      unsigned direction = edges_of_verts_directions[j];
      if (!(col_codes[edge] & (1<<direction)))
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
  *candidates_out = candidates;
  *gen_vert_of_verts_out = gen_vert_of_verts;
  *col_qual_of_verts_out = col_qual_of_verts;
}
