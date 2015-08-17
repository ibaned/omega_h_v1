#include "coarsen_by_size.h"
#include "collapse_codes.h"
#include "measure_edges.h"
#include "ints.h"
#include "check_collapse_class.h"
#include "coarsen_qualities.h"
#include "collapses_to_verts.h"
#include "indset.h"
#include "collapses_to_elements.h"
#include "coarsen_topology.h"
#include "concat.h"
#include "tables.h"
#include "element_qualities.h"
#include <stdlib.h>

unsigned coarsen_by_size(
    struct mesh** p_m,
    double (*size_function)(double const x[]),
    double quality_floor)
{
  struct mesh* m = *p_m;
  unsigned nedges = mesh_count(m, 1);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  unsigned nverts = mesh_count(m, 0);
  double* sizes = malloc(sizeof(double) * nverts);
  for (unsigned i = 0; i < nverts; ++i)
    sizes[i] = size_function(coords + i * 3);
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, sizes);
  free(sizes);
  unsigned* col_codes = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    if (edge_sizes[i] < (1.0 / 2.0))
      col_codes[i] = COLLAPSE_BOTH;
    else
      col_codes[i] = DONT_COLLAPSE;
  }
  free(edge_sizes);
  if (ints_max(col_codes, nedges) == DONT_COLLAPSE) {
    /* early return #1: no edges are small */
    free(col_codes);
    return 0;
  }
  unsigned elem_dim = mesh_dim(m);
  unsigned const* class_dim = mesh_find_nodal_label(m, "class_dim")->data;
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
  double* quals_of_edges = coarsen_qualities(elem_dim, nedges, col_codes,
      verts_of_elems, verts_of_edges, elems_of_verts_offsets,
      elems_of_verts, elems_of_verts_directions, coords, quality_floor);
  if (ints_max(col_codes, nedges) == DONT_COLLAPSE) {
    /* early return #2: all small edges failed their classif/quality checks */
    free(col_codes);
    free(quals_of_edges);
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
  free(quals_of_edges);
  free(col_codes);
  unsigned* indset = find_indset(nverts, verts_of_verts_offsets, verts_of_verts,
      candidates, qual_of_verts);
  free(candidates);
  free(qual_of_verts);
  unsigned* gen_offset_of_verts = ints_exscan(indset, nverts);
  free(indset);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned* gen_offset_of_elems;
  unsigned* gen_vert_of_elems;
  unsigned* gen_direction_of_elems;
  unsigned* offset_of_same_elems;
  collapses_to_elements(elem_dim, nelems, verts_of_elems, gen_offset_of_verts,
      gen_vert_of_verts, &gen_offset_of_elems, &gen_vert_of_elems,
      &gen_direction_of_elems, &offset_of_same_elems);
  free(gen_vert_of_verts);
  unsigned* offset_of_same_verts = ints_negate_offsets(
      gen_offset_of_verts, nverts);
  unsigned nverts_out = offset_of_same_verts[nverts];
  free(gen_offset_of_verts);
  double* coords_out = doubles_subset(nverts, 3, coords, offset_of_same_verts);
  unsigned* class_dim_out = ints_subset(nverts, 1, class_dim,
      offset_of_same_verts);
  unsigned ngen_elems;
  unsigned* verts_of_gen_elems;
  coarsen_topology(elem_dim, nelems, verts_of_elems, gen_offset_of_elems,
      gen_vert_of_elems, gen_direction_of_elems, &ngen_elems,
      &verts_of_gen_elems);
  free(gen_offset_of_elems);
  free(gen_vert_of_elems);
  free(gen_direction_of_elems);
  unsigned nelems_out;
  unsigned* verts_of_elems_out;
  concat_verts_of_elems(elem_dim, nelems, ngen_elems, verts_of_elems,
      offset_of_same_elems, verts_of_gen_elems,
      &nelems_out, &verts_of_elems_out);
  free(offset_of_same_elems);
  free(verts_of_gen_elems);
  /* remap element vertices to account for vertex removal */
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  for (unsigned i = 0; i < nelems_out * verts_per_elem; ++i)
    verts_of_elems_out[i] = offset_of_same_verts[verts_of_elems_out[i]];
  free(offset_of_same_verts);
  struct mesh* m_out = new_mesh(elem_dim);
  mesh_set_ents(m_out, 0, nverts_out, 0);
  mesh_set_ents(m_out, elem_dim, nelems_out, verts_of_elems_out);
  mesh_add_nodal_field(m_out, "coordinates", 3, coords_out);
  mesh_add_nodal_label(m_out, "class_dim", class_dim_out);
  free_mesh(m);
  *p_m = m_out;
  return 1;
}
