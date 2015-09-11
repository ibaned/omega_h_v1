#include "quality.h"
#include "swap_qualities.h"
#include "swap_topology.h"
#include "mark.h"
#include "ints.h"
#include "mesh.h"
#include "concat.h"
#include <stdlib.h>
#include <assert.h>

unsigned swap_common(
    struct mesh** p_m,
    unsigned* candidates,
    double good_qual,
    double valid_qual,
    unsigned require_better)
{
  struct mesh* m = *p_m;
  assert(mesh_dim(m) == 3);
  unsigned nedges = mesh_count(m, 1);
  if (!ints_max(candidates, nedges))
    return 0;
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* class_dim = mesh_find_nodal_label(m, "class_dim")->data;
  unmark_boundary(3, 1, nedges, verts_of_edges, class_dim, candidates);
  if (!ints_max(candidates, nedges))
    return 0;
  unsigned const* elems_of_edges_offsets =
    mesh_ask_up(m, 1, 3)->offsets;
  unsigned const* elems_of_edges =
    mesh_ask_up(m, 1, 3)->adj;
  unsigned const* elems_of_edges_directions =
    mesh_ask_up(m, 1, 3)->directions;
  unsigned const* verts_of_elems = mesh_ask_down(m, 3, 0);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  double* elem_quals = 0;
  if (require_better)
    elem_quals = mesh_qualities(m);
  double* edge_quals;
  unsigned* edge_codes;
  unsigned* gen_elems_per_edge;
  swap_qualities(nedges, candidates,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions,
      verts_of_edges, verts_of_elems, coords,
      good_qual, valid_qual, elem_quals, require_better,
      &edge_quals, &edge_codes, &gen_elems_per_edge);
  free(elem_quals);
  if (!ints_max(candidates, nedges)) {
    free(edge_quals);
    free(edge_codes);
    free(gen_elems_per_edge);
    return 0;
  }
  unsigned* gen_offset_of_edges = ints_exscan(gen_elems_per_edge, nedges);
  free(gen_elems_per_edge);
  unsigned* verts_of_gen_elems = swap_topology(nedges, candidates,
      gen_offset_of_edges, edge_codes,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions,
      verts_of_edges, verts_of_elems);
  unsigned ngen_elems = gen_offset_of_edges[nedges];
  unsigned* old_elems = mesh_mark_up(m, 1, 3, candidates);
  unsigned nelems = mesh_count(m, 3);
  unsigned* same_elems = ints_negate(old_elems, nelems);
  free(old_elems);
  unsigned* same_elem_offsets = ints_exscan(same_elems, nelems);
  free(same_elems);
  unsigned nelems_out;
  unsigned* verts_of_elems_out;
  concat_verts_of_elems(3, nelems, ngen_elems, verts_of_elems,
      same_elem_offsets, verts_of_gen_elems,
      &nelems_out, &verts_of_elems_out);
  struct mesh* m_out = new_mesh(3);
  /* todo... */
  (void) m_out;
  return 1;
}
