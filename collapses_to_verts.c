#include "collapses_to_verts.h"

#include "collapse_codes.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"

void valid_collapses_to_verts(struct mesh* m)
{
  unsigned const* col_codes = mesh_find_tag(m, 1, "col_codes");
  double const* col_quals_of_edges =
    mesh_find_tag(m, 1, "quals_of_edges")->d.f64;
  unsigned nverts = mesh_count(m, 0);
  unsigned const* edges_of_verts_offsets =
    mesh_ask_up(m, 0, 1)->offsets;
  unsigned const* edges_of_verts =
    mesh_ask_up(m, 0, 1)->adj;
  unsigned const* edges_of_verts_directions =
    mesh_ask_up(m, 0, 1)->directions;
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
  if (mesh_is_parallel(m)) {
    mesh_conform_uints(m, 0, 1, &candidates);
    mesh_conform_doubles(m, 0, 1, &col_qual_of_verts);
  }
  mesh_add_tag(m, 0, TAG_U32, "candidates", 1, candidates);
  mesh_add_tag(m, 0, TAG_U32, "col_qual", 1, col_qual_of_verts);
}
