#include "collapses_to_verts.h"
#include "tables.h"
#include "ints.h"
#include <stdlib.h>

void collapses_to_verts(
    unsigned nverts,
    unsigned const* verts_of_edges,
    unsigned const* edges_of_verts_offsets,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_directions,
    unsigned const* col_codes,
    double const* col_quals_of_edges,
    unsigned** gen_offset_of_verts_out,
    unsigned** gen_vert_of_verts_out,
    double** col_qual_of_verts_out)
{
  unsigned* col_verts = malloc(sizeof(unsigned) * nverts);
  unsigned* gen_vert_of_verts = malloc(sizeof(unsigned) * nverts);
  double* col_qual_of_verts = malloc(sizeof(double) * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_use = edges_of_verts_offsets[i];
    unsigned end_use = edges_of_verts_offsets[i + 1];
    double maxq;
    unsigned gen_vert;
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
    col_verts[i] = is_collapsing;
    if (is_collapsing) {
      gen_vert_of_verts[i] = gen_vert;
      col_qual_of_verts[i] = maxq;
    }
  }
  *gen_offset_of_verts_out = ints_exscan(col_verts, nverts);
  *gen_vert_of_verts_out = gen_vert_of_verts;
  *col_qual_of_verts_out = col_qual_of_verts;
}
