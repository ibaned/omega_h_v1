#include "coarsen_reduced.h"
#include "collapse_codes.h"
#include "up_from_down.h"
#include "star.h"
#include "bridge_graph.h"
#include "measure_edges.h"
#include "ints.h"
#include "reflect_down.h"
#include "collapse_classif.h"
#include "coarsen_qualities.h"
#include "collapses_to_verts.h"
#include "indset.h"
#include "collapses_to_elements.h"
#include "coarsen_topology.h"
#include "concat.h"
#include <stdlib.h>

int coarsen_reduced(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    unsigned** p_class_dim,
    double (*size_function)(double const x[]))
{
  unsigned nelems = *p_nelems;
  unsigned nverts = *p_nverts;
  unsigned const* verts_of_elems = *p_verts_of_elems;
  double const* coords = *p_coords;
  unsigned const* class_dim = *p_class_dim;
  unsigned* elems_of_verts_offsets;
  unsigned* elems_of_verts;
  unsigned* elems_of_verts_directions;
  up_from_down(elem_dim, 0, nelems, nverts, verts_of_elems,
      &elems_of_verts_offsets, &elems_of_verts, &elems_of_verts_directions);
  unsigned* verts_of_verts_offsets;
  unsigned* verts_of_verts;
  get_star(0, elem_dim, nverts, elems_of_verts_offsets, elems_of_verts,
      verts_of_elems, &verts_of_verts_offsets, &verts_of_verts);
  unsigned nedges;
  unsigned* verts_of_edges;
  bridge_graph(nverts, verts_of_verts_offsets, verts_of_verts,
      &nedges, &verts_of_edges);
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
    free(verts_of_verts_offsets);
    free(verts_of_verts);
    free(verts_of_edges);
    free(col_codes);
    return 0;
  }
  unsigned* edges_of_verts_offsets;
  unsigned* edges_of_verts;
  unsigned* edges_of_verts_directions;
  up_from_down(1, 0, nedges, nverts, verts_of_edges,
      &edges_of_verts_offsets, &edges_of_verts, &edges_of_verts_directions);
  unsigned* edges_of_elems = reflect_down(elem_dim, 1, nelems, nedges,
      verts_of_elems, edges_of_verts_offsets, edges_of_verts);
  unsigned* elems_of_edges_offsets;
  unsigned* elems_of_edges;
  unsigned* elems_of_edges_directions;
  up_from_down(elem_dim, 1, nelems, nedges, edges_of_elems,
      &elems_of_edges_offsets, &elems_of_edges, &elems_of_edges_directions);
  free(edges_of_elems);
  check_collapse_classif(elem_dim, nedges, col_codes, class_dim,
      verts_of_elems, verts_of_edges, verts_of_verts_offsets, verts_of_verts,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions);
  free(elems_of_edges_offsets);
  free(elems_of_edges);
  free(elems_of_edges_directions);
  double* quals_of_edges = coarsen_qualities(elem_dim, nedges, col_codes,
      verts_of_elems, verts_of_edges, elems_of_verts_offsets,
      elems_of_verts, elems_of_verts_directions, coords);
  free(elems_of_verts_offsets);
  free(elems_of_verts);
  free(elems_of_verts_directions);
  if (ints_max(col_codes, nedges) == DONT_COLLAPSE) {
    /* early return #2: all small edges failed their classif/quality checks */
    free(verts_of_verts_offsets);
    free(verts_of_verts);
    free(verts_of_edges);
    free(col_codes);
    free(edges_of_verts_offsets);
    free(edges_of_verts);
    free(edges_of_verts_directions);
    free(quals_of_edges);
    return 0;
  }
  /* from this point forward, some edges will definitely collapse */
  unsigned* candidates;
  unsigned* gen_vert_of_verts;
  double* qual_of_verts;
  collapses_to_verts(nverts, verts_of_edges, edges_of_verts_offsets,
      edges_of_verts, edges_of_verts_directions, col_codes, quals_of_edges,
      &candidates, &gen_vert_of_verts, &qual_of_verts);
  free(edges_of_verts_offsets);
  free(edges_of_verts);
  free(edges_of_verts_directions);
  free(verts_of_edges);
  free(quals_of_edges);
  free(col_codes);
  unsigned* indset = find_indset(nverts, verts_of_verts_offsets, verts_of_verts,
      candidates, qual_of_verts);
  free(verts_of_verts_offsets);
  free(verts_of_verts);
  free(candidates);
  free(qual_of_verts);
  unsigned* gen_offset_of_verts = ints_exscan(indset, nverts);
  free(indset);
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
  free(offset_of_same_verts);
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
  *p_nelems = nelems_out;
  *p_nverts = nverts_out;
  free(*p_verts_of_elems);
  *p_verts_of_elems = verts_of_elems_out;
  free(*p_coords);
  *p_coords = coords_out;
  free(*p_class_dim);
  *p_class_dim = class_dim_out;
  return 1;
}
