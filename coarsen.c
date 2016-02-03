#include "coarsen.h"

#include "coarsen_common.h"
#include "collapse_codes.h"
#include "ghost_mesh.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "size.h"

LOOP_KERNEL(coarsen_size_code,
    double const* edge_sizes,
    double size_ratio_floor,
    unsigned* col_codes)
  if (edge_sizes[i] < size_ratio_floor)
    col_codes[i] = COLLAPSE_BOTH;
  else
    col_codes[i] = DONT_COLLAPSE;
}

unsigned coarsen_by_size(
    struct mesh** p_m,
    double quality_floor,
    double size_ratio_floor)
{
  struct mesh* m = *p_m;
  unsigned nedges = mesh_count(m, 1);
  double* edge_sizes = mesh_measure_edges_for_adapt(m);
  unsigned* col_codes = LOOP_MALLOC(unsigned, nedges);
  LOOP_EXEC(coarsen_size_code, nedges,
      edge_sizes, size_ratio_floor, col_codes);
  loop_free(edge_sizes);
  mesh_add_tag(m, 1, TAG_U32, "col_codes", 1, col_codes);
  return coarsen_common(p_m, quality_floor, 0);
}

LOOP_KERNEL(coarsen_sliver_code,
    unsigned const* verts_of_edges,
    unsigned const* marked_verts,
    unsigned* col_codes)
  unsigned const* verts_of_edge = verts_of_edges + i * 2;
  col_codes[i] = 0;
  for (unsigned j = 0; j < 2; ++j) {
    unsigned vert = verts_of_edge[j];
    if (marked_verts[vert])
      col_codes[i] = do_collapse(col_codes[i], j);
  }
}

unsigned coarsen_slivers(
    struct mesh** p_m,
    double quality_floor,
    unsigned nlayers)
{
  if (mesh_is_parallel(*p_m))
    mesh_ensure_ghosting(*p_m, 1);
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned* slivers = mesh_mark_slivers(m, quality_floor, nlayers);
  unsigned* marked_verts = mesh_mark_down(m, elem_dim, 0, slivers);
  loop_free(slivers);
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned* col_codes = LOOP_MALLOC(unsigned, nedges);
  LOOP_EXEC(coarsen_sliver_code, nedges,
      verts_of_edges, marked_verts, col_codes);
  loop_free(marked_verts);
  mesh_add_tag(m, 1, TAG_U32, "col_codes", 1, col_codes);
  return coarsen_common(p_m, 0.0, 1);
}
