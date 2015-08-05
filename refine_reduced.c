#include "refine_reduced.h"
#include "derive_edges.h"
#include "measure_edges.h"
#include "refine_nodal.h"
#include "up_from_down.h"
#include "star.h"
#include "reflect_down.h"
#include "indset.h"
#include "ints.h"
#include "splits_to_elements.h"
#include "refine_topology.h"
#include "concat.h"
#include <stdlib.h>

static void derive_adjacencies(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    unsigned* nedges_out,
    unsigned** verts_of_edges_out,
    unsigned** edges_of_verts_offsets_out,
    unsigned** edges_of_verts_out,
    unsigned** edges_of_elems_out,
    unsigned** elems_of_edges_offsets_out,
    unsigned** elems_of_edges_out)
{
  derive_edges(elem_dim, nelems, nverts,
      verts_of_elems,
      nedges_out, verts_of_edges_out);
  unsigned nedges = *nedges_out;
  unsigned const* verts_of_edges = *verts_of_edges_out;
  up_from_down(1, 0, nedges, nverts,
      verts_of_edges,
      edges_of_verts_offsets_out, edges_of_verts_out, 0);
  unsigned const* edges_of_verts_offsets = *edges_of_verts_offsets_out;
  unsigned const* edges_of_verts = *edges_of_verts_out;
  *edges_of_elems_out = reflect_down(elem_dim, 1, nelems, nedges,
      verts_of_elems,
      edges_of_verts_offsets, edges_of_verts);
  unsigned const* edges_of_elems = *edges_of_elems_out;
  up_from_down(elem_dim, 1, nelems, nedges,
      edges_of_elems,
      elems_of_edges_offsets_out, elems_of_edges_out, 0);
}

static unsigned* choose_edges(
    unsigned elem_dim,
    unsigned nedges,
    unsigned nverts,
    unsigned const* verts_of_edges,
    unsigned const* edges_of_elems,
    unsigned const* elems_of_edges_offsets,
    unsigned const* elems_of_edges,
    double const* coords,
    double (*size_function)(double const x[]))
{
  double* sizes = malloc(sizeof(double) * nverts);
  for (unsigned i = 0; i < nverts; ++i)
    sizes[i] = size_function(coords + i * 3);
  double* edge_sizes = measure_edges(nedges, verts_of_edges, coords, sizes);
  free(sizes);
  unsigned* candidates = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i)
    candidates[i] = edge_sizes[i] > 1.0;
  unsigned* edges_of_edges_offsets;
  unsigned* edges_of_edges;
  get_star(1, elem_dim, nedges,
      elems_of_edges_offsets, elems_of_edges,
      edges_of_elems,
      &edges_of_edges_offsets, &edges_of_edges);
  unsigned* indset = find_indset(nedges,
      edges_of_edges_offsets, edges_of_edges,
      candidates,
      edge_sizes);
  free(edges_of_edges_offsets);
  free(edges_of_edges);
  free(candidates);
  free(edge_sizes);
  return indset;
}

struct reduced_mesh refine_reduced(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords,
    double (*size_function)(double const x[]))
{
  unsigned nedges;
  unsigned* verts_of_edges;
  unsigned* edges_of_verts_offsets;
  unsigned* edges_of_verts;
  unsigned* edges_of_elems;
  unsigned* elems_of_edges_offsets;
  unsigned* elems_of_edges;
  derive_adjacencies(elem_dim, nelems, nverts, verts_of_elems,
      &nedges,
      &verts_of_edges,
      &edges_of_verts_offsets, &edges_of_verts,
      &edges_of_elems,
      &elems_of_edges_offsets, &elems_of_edges);
  unsigned* indset = choose_edges(elem_dim, nedges, nverts, verts_of_edges,
      edges_of_elems,
      elems_of_edges_offsets, elems_of_edges,
      coords,
      size_function);
  unsigned* gen_offset_of_edges = ints_exscan(indset, nedges);
  free(indset);
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
  unsigned ngen_elems;
  unsigned* verts_of_gen_elems;
  refine_topology(elem_dim, 1, elem_dim, nelems, verts_of_elems,
      gen_offset_of_elems, gen_vert_of_elems, gen_direction_of_elems,
      &ngen_elems, &verts_of_gen_elems);
  free(gen_offset_of_elems);
  free(gen_vert_of_elems);
  free(gen_direction_of_elems);
  double* gen_coords = refine_nodal(1, nedges, verts_of_edges,
      gen_offset_of_edges, 3, coords);
  unsigned concat_sizes[2] = {nverts, nedges};
  double const* concat_coords[2] = {coords, gen_coords};
  double* new_coords = concat_doubles(2, 3, concat_sizes, concat_coords);
  free(gen_coords);
  (void) new_coords;
  struct reduced_mesh rm;
  return rm;
}
