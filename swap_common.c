#include "swap_common.h"
#include "quality.h"
#include "swap_qualities.h"
#include "swap_topology.h"
#include "mark.h"
#include "ints.h"
#include "doubles.h"
#include "mesh.h"
#include "concat.h"
#include "indset.h"
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
  unsigned elem_dim = mesh_dim(m);
  assert(elem_dim == 3);
  unsigned nedges = mesh_count(m, 1);
  if (!ints_max(candidates, nedges))
    return 0;
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* class_dim = mesh_find_nodal_label(m, "class_dim")->data;
  unmark_boundary(elem_dim, 1, nedges, verts_of_edges, class_dim, candidates);
  if (!ints_max(candidates, nedges))
    return 0;
  unsigned const* elems_of_edges_offsets =
    mesh_ask_up(m, 1, elem_dim)->offsets;
  unsigned const* elems_of_edges =
    mesh_ask_up(m, 1, elem_dim)->adj;
  unsigned const* elems_of_edges_directions =
    mesh_ask_up(m, 1, elem_dim)->directions;
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
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
  unsigned const* edges_of_edges_offsets = mesh_ask_star(m, 1, elem_dim)->offsets;
  unsigned const* edges_of_edges = mesh_ask_star(m, 1, elem_dim)->adj;
  unsigned* indset = find_indset(nedges, edges_of_edges_offsets, edges_of_edges,
      candidates, edge_quals);
  for (unsigned i = 0; i < nedges; ++i)
    if (!indset[i])
      gen_elems_per_edge[i] = 0;
  unsigned* gen_offset_of_edges = ints_exscan(gen_elems_per_edge, nedges);
  free(gen_elems_per_edge);
  unsigned* verts_of_gen_elems = swap_topology(nedges, indset,
      gen_offset_of_edges, edge_codes,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions,
      verts_of_edges, verts_of_elems);
  unsigned ngen_elems = gen_offset_of_edges[nedges];
  unsigned* old_elems = mesh_mark_up(m, 1, elem_dim, indset);
  free(indset);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned* same_elems = ints_negate(old_elems, nelems);
  free(old_elems);
  unsigned* same_elem_offsets = ints_exscan(same_elems, nelems);
  free(same_elems);
  unsigned nelems_out;
  unsigned* verts_of_elems_out;
  concat_verts_of_elems(elem_dim, nelems, ngen_elems, verts_of_elems,
      same_elem_offsets, verts_of_gen_elems,
      &nelems_out, &verts_of_elems_out);
  unsigned nverts = mesh_count(m, 0);
  struct mesh* m_out = new_mesh(elem_dim);
  mesh_set_ents(m_out, 0, nverts, 0);
  mesh_set_ents(m_out, elem_dim, nelems_out, verts_of_elems_out);
  for (unsigned i = 0; i < mesh_count_nodal_fields(m); ++i) {
    struct const_field* f = mesh_get_nodal_field(m, i);
    double* data = doubles_copy(f->data, nverts * f->ncomps);
    mesh_add_nodal_field(m_out, f->name, f->ncomps, data);
  }
  for (unsigned i = 0; i < mesh_count_nodal_labels(m); ++i) {
    struct const_label* l = mesh_get_nodal_label(m, i);
    unsigned* data = ints_copy(l->data, nverts);
    mesh_add_nodal_label(m_out, l->name, data);
  }
  free_mesh(m);
  *p_m = m_out;
  return 1;
}
