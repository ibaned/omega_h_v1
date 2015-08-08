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
  unsigned something_to_do = ints_max(col_codes, nedges);
  if (!something_to_do) {
    free(verts_of_edges);
    free(verts_of_verts_offsets);
    free(verts_of_verts);
    return 0;
  }
  unsigned* edges_of_verts_offsets;
  unsigned* edges_of_verts;
  unsigned* edges_of_verts_directions;
  up_from_down(1, 0, nedges, nverts, verts_of_edges,
      &edges_of_verts_offsets, &edges_of_verts, &edges_of_verts_directions);
  unsigned* edges_of_elems = reflect_down(elem_dim, 1, nelems, nedges,
      verts_of_elems, edges_of_verts_offsets, edges_of_verts);
  free(edges_of_verts_offsets);
  free(edges_of_verts);
  unsigned* elems_of_edges_offsets;
  unsigned* elems_of_edges;
  unsigned* elems_of_edges_directions;
  up_from_down(elem_dim, 1, nelems, nedges, edges_of_elems,
      &elems_of_edges_offsets, &elems_of_edges, &elems_of_edges_directions);
  check_collapse_classif(elem_dim, nedges, col_codes, class_dim,
      verts_of_elems, verts_of_edges, verts_of_verts_offsets, verts_of_verts,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions);
  double* quals_of_edges = coarsen_qualities(elem_dim, nedges, col_codes,
      verts_of_elems, verts_of_edges, elems_of_verts_offsets,
      elems_of_verts, elems_of_verts_directions, coords);
  unsigned* gen_offset_of_verts;
  unsigned* gen_vert_of_verts;
  double* qual_of_verts;
  collapses_to_verts(nverts, verts_of_edges, edges_of_verts_offsets,
      edges_of_verts, edges_of_verts_directions, col_codes, quals_of_edges,
      &gen_offset_of_verts, &gen_vert_of_verts, &qual_of_verts);
  free(quals_of_edges);
  free(col_codes);
  unsigned* candidates = ints_unscan(gen_offset_of_verts, nverts);
  unsigned* indset = find_indset(nverts, verts_of_verts_offsets, verts_of_verts,
      candidates, qual_of_verts);
  free(candidates);
  free(qual_of_verts);
  free(gen_offset_of_verts);
  gen_offset_of_verts = ints_exscan(indset, nverts);
  unsigned* gen_offset_of_elems;
  unsigned* gen_vert_of_elems;
  unsigned* gen_direction_of_elems;
  collapses_to_elements(elem_dim, nelems, verts_of_elems, gen_offset_of_verts,
      gen_vert_of_verts, &gen_offset_of_elems, &gen_vert_of_elems,
      &gen_direction_of_elems);
  unsigned ngen_elems;
  unsigned* verts_of_gen_elems;
  coarsen_topology(elem_dim, nelems, verts_of_elems, gen_offset_of_elems,
      gen_vert_of_elems, gen_direction_of_elems, &ngen_elems,
      &verts_of_gen_elems);
  return 1;
}
