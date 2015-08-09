#include "refine_reduced.h"
#include "measure_edges.h"
#include "refine_nodal.h"
#include "up_from_down.h"
#include "star.h"
#include "reflect_down.h"
#include "indset.h"
#include "ints.h"
#include "splits_to_elements.h"
#include "refine_topology.h"
#include "refine_qualities.h"
#include "concat.h"
#include "tables.h"
#include "bridge_graph.h"
#include <stdlib.h>

int refine_reduced(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    double (*size_function)(double const x[]))
{
  unsigned nelems = *p_nelems;
  unsigned nverts = *p_nverts;
  unsigned const* verts_of_elems = *p_verts_of_elems;
  double const* coords = *p_coords;
  unsigned* elems_of_verts_offsets;
  unsigned* elems_of_verts;
  up_from_down(elem_dim, 0, nelems, nverts, verts_of_elems,
      &elems_of_verts_offsets, &elems_of_verts, 0);
  unsigned* verts_of_verts_offsets;
  unsigned* verts_of_verts;
  get_star(0, elem_dim, nverts, elems_of_verts_offsets, elems_of_verts,
      verts_of_elems, &verts_of_verts_offsets, &verts_of_verts);
  free(elems_of_verts_offsets);
  free(elems_of_verts);
  unsigned nedges;
  unsigned* verts_of_edges;
  bridge_graph(nverts, verts_of_verts_offsets, verts_of_verts,
      &nedges, &verts_of_edges);
  free(verts_of_verts_offsets);
  free(verts_of_verts);
  double* sizes = malloc(sizeof(double) * nverts);
  for (unsigned i = 0; i < nverts; ++i)
    sizes[i] = size_function(coords + i * 3);
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, sizes);
  free(sizes);
  unsigned* candidates = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i)
    candidates[i] = edge_sizes[i] > 1.0;
  free(edge_sizes);
  unsigned something_to_do = ints_max(candidates, nedges);
  if (!something_to_do) {
    free(verts_of_edges);
    free(candidates);
    return 0;
  }
  unsigned* edges_of_verts_offsets;
  unsigned* edges_of_verts;
  up_from_down(1, 0, nedges, nverts, verts_of_edges,
      &edges_of_verts_offsets, &edges_of_verts, 0);
  unsigned* edges_of_elems = reflect_down(elem_dim, 1, nelems, nedges,
      verts_of_elems, edges_of_verts_offsets, edges_of_verts);
  free(edges_of_verts_offsets);
  free(edges_of_verts);
  unsigned* elems_of_edges_offsets;
  unsigned* elems_of_edges;
  unsigned* elems_of_edges_directions;
  up_from_down(elem_dim, 1, nelems, nedges, edges_of_elems,
      &elems_of_edges_offsets, &elems_of_edges, &elems_of_edges_directions);
  double* edge_quals = refine_qualities(elem_dim, 1, nedges, verts_of_edges,
      verts_of_elems, elems_of_edges_offsets, elems_of_edges,
      elems_of_edges_directions, candidates, coords);
  free(elems_of_edges_directions);
  unsigned* edges_of_edges_offsets;
  unsigned* edges_of_edges;
  get_star(1, elem_dim, nedges, elems_of_edges_offsets, elems_of_edges,
      edges_of_elems, &edges_of_edges_offsets, &edges_of_edges);
  unsigned* indset = find_indset(nedges, edges_of_edges_offsets, edges_of_edges,
      candidates, edge_quals);
  free(elems_of_edges_offsets);
  free(elems_of_edges);
  free(edges_of_edges_offsets);
  free(edges_of_edges);
  free(candidates);
  free(edge_quals);
  unsigned* gen_offset_of_edges = ints_exscan(indset, nedges);
  free(indset);
  unsigned nsplit_edges = gen_offset_of_edges[nedges];
  unsigned* gen_vert_of_edges = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i)
    if (gen_offset_of_edges[i] != gen_offset_of_edges[i + 1])
      gen_vert_of_edges[i] = nverts + gen_offset_of_edges[i];
  unsigned* gen_offset_of_elems;
  unsigned* gen_direction_of_elems;
  unsigned* gen_vert_of_elems;
  project_splits_to_elements(elem_dim, 1, nelems,
      edges_of_elems, gen_offset_of_edges, gen_vert_of_edges,
      &gen_offset_of_elems, &gen_direction_of_elems, &gen_vert_of_elems);
  free(edges_of_elems);
  free(gen_vert_of_edges);
  unsigned ngen_elems;
  unsigned* verts_of_gen_elems;
  refine_topology(elem_dim, 1, elem_dim, nelems, verts_of_elems,
      gen_offset_of_elems, gen_vert_of_elems, gen_direction_of_elems,
      &ngen_elems, &verts_of_gen_elems);
  free(gen_vert_of_elems);
  free(gen_direction_of_elems);
  double* gen_coords = refine_nodal(1, nedges, verts_of_edges,
      gen_offset_of_edges, 3, coords);
  free(verts_of_edges);
  free(gen_offset_of_edges);
  unsigned concat_sizes[2] = {nverts, nsplit_edges};
  unsigned nverts_out = nverts + nsplit_edges;
  double const* concat_coords[2] = {coords, gen_coords};
  double* coords_out = concat_doubles(2, 3, concat_sizes, concat_coords);
  free(gen_coords);
  unsigned* offset_of_same_elems = ints_negate_offsets(
      gen_offset_of_elems, nelems);
  free(gen_offset_of_elems);
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
  return 1;
}
