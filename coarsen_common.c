#include "coarsen_common.h"
#include "loop.h"                 // for free
#include "check_collapse_class.h"   // for check_collapse_class
#include "coarsen_qualities.h"      // for coarsen_qualities
#include "coarsen_topology.h"       // for coarsen_topology
#include "collapse_codes.h"         // for ::DONT_COLLAPSE
#include "collapses_to_elements.h"  // for collapses_to_elements
#include "collapses_to_verts.h"     // for collapses_to_verts
#include "concat.h"                 // for concat_verts_of_elems
#include "field.h"                  // for const_field
#include "graph.h"                  // for const_graph
#include "indset.h"                 // for find_indset
#include "ints.h"                   // for ints_max, ints_exscan, ints_negat...
#include "label.h"                  // for const_label
#include "mesh.h"                   // for mesh_ask_up, const_up, mesh_count
#include "quality.h"                // for element_qualities
#include "subset.h"                 // for doubles_subset, ints_subset
#include "tables.h"                 // for the_down_degrees

unsigned coarsen_common(
    struct mesh** p_m,
    unsigned* col_codes,
    double quality_floor,
    unsigned require_better)
{
  struct mesh* m = *p_m;
  unsigned nedges = mesh_count(m, 1);
  if (ints_max(col_codes, nedges) == DONT_COLLAPSE)
    return 0;
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  double const* coords = mesh_find_field(m, 0, "coordinates")->data;
  unsigned nverts = mesh_count(m, 0);
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* class_dim = mesh_find_label(m, 0, "class_dim")->data;
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned const* verts_of_verts_offsets = mesh_ask_star(m, 0, elem_dim)->offsets;
  unsigned const* verts_of_verts = mesh_ask_star(m, 0, elem_dim)->adj;
  unsigned const* elems_of_edges_offsets = mesh_ask_up(m, 1, elem_dim)->offsets;
  unsigned const* elems_of_edges = mesh_ask_up(m, 1, elem_dim)->adj;
  unsigned const* elems_of_edges_directions =
    mesh_ask_up(m, 1, elem_dim)->directions;
  check_collapse_class(elem_dim, nedges, col_codes, class_dim,
      verts_of_elems, verts_of_edges, verts_of_verts_offsets, verts_of_verts,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions);
  unsigned const* elems_of_verts_offsets = mesh_ask_up(m, 0, elem_dim)->offsets;
  unsigned const* elems_of_verts = mesh_ask_up(m, 0, elem_dim)->adj;
  unsigned const* elems_of_verts_directions = mesh_ask_up(m, 0, elem_dim)->directions;
  double* elem_quals = 0;
  if (require_better)
    elem_quals = element_qualities(elem_dim, nelems, verts_of_elems, coords);
  double* quals_of_edges = coarsen_qualities(elem_dim, nedges, col_codes,
      verts_of_elems, verts_of_edges, elems_of_verts_offsets,
      elems_of_verts, elems_of_verts_directions, coords, quality_floor,
      elem_quals, require_better);
  loop_free(elem_quals);
  if (ints_max(col_codes, nedges) == DONT_COLLAPSE) {
    /* early return #2: all small edges failed their classif/quality checks */
    loop_free(quals_of_edges);
    return 0;
  }
  /* from this point forward, some edges will definitely collapse */
  unsigned* candidates;
  unsigned* gen_vert_of_verts;
  double* qual_of_verts;
  unsigned const* edges_of_verts_offsets = mesh_ask_up(m, 0, 1)->offsets;
  unsigned const* edges_of_verts = mesh_ask_up(m, 0, 1)->adj;
  unsigned const* edges_of_verts_directions = mesh_ask_up(m, 0, 1)->directions;
  collapses_to_verts(nverts, verts_of_edges, edges_of_verts_offsets,
      edges_of_verts, edges_of_verts_directions, col_codes, quals_of_edges,
      &candidates, &gen_vert_of_verts, &qual_of_verts);
  loop_free(quals_of_edges);
  unsigned* indset = find_indset(nverts, verts_of_verts_offsets, verts_of_verts,
      candidates, qual_of_verts);
  loop_free(candidates);
  loop_free(qual_of_verts);
  unsigned* gen_offset_of_verts = ints_exscan(indset, nverts);
  loop_free(indset);
  unsigned* gen_offset_of_elems;
  unsigned* gen_vert_of_elems;
  unsigned* gen_direction_of_elems;
  unsigned* offset_of_same_elems;
  collapses_to_elements(elem_dim, nelems, verts_of_elems, gen_offset_of_verts,
      gen_vert_of_verts, &gen_offset_of_elems, &gen_vert_of_elems,
      &gen_direction_of_elems, &offset_of_same_elems);
  loop_free(gen_vert_of_verts);
  unsigned* offset_of_same_verts = ints_negate_offsets(
      gen_offset_of_verts, nverts);
  unsigned nverts_out = offset_of_same_verts[nverts];
  loop_free(gen_offset_of_verts);
  unsigned ngen_elems;
  unsigned* verts_of_gen_elems;
  coarsen_topology(elem_dim, nelems, verts_of_elems, gen_offset_of_elems,
      gen_vert_of_elems, gen_direction_of_elems, &ngen_elems,
      &verts_of_gen_elems);
  loop_free(gen_offset_of_elems);
  loop_free(gen_vert_of_elems);
  loop_free(gen_direction_of_elems);
  unsigned nelems_out;
  unsigned* verts_of_elems_out;
  concat_verts_of_elems(elem_dim, nelems, ngen_elems, verts_of_elems,
      offset_of_same_elems, verts_of_gen_elems,
      &nelems_out, &verts_of_elems_out);
  loop_free(offset_of_same_elems);
  loop_free(verts_of_gen_elems);
  /* remap element vertices to account for vertex removal */
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  for (unsigned i = 0; i < nelems_out * verts_per_elem; ++i)
    verts_of_elems_out[i] = offset_of_same_verts[verts_of_elems_out[i]];
  struct mesh* m_out = new_mesh(elem_dim);
  mesh_set_ents(m_out, 0, nverts_out, 0);
  mesh_set_ents(m_out, elem_dim, nelems_out, verts_of_elems_out);
  for (unsigned i = 0; i < mesh_count_fields(m, 0); ++i) {
    struct const_field* f = mesh_get_field(m, 0, i);
    double* vals_out = doubles_subset(nverts, f->ncomps, f->data,
        offset_of_same_verts);
    mesh_add_field(m_out, 0, f->name, f->ncomps, vals_out);
  }
  unsigned* class_dim_out = ints_subset(nverts, 1, class_dim,
      offset_of_same_verts);
  mesh_add_label(m_out, 0, "class_dim", class_dim_out);
  loop_free(offset_of_same_verts);
  free_mesh(m);
  *p_m = m_out;
  return 1;
}
