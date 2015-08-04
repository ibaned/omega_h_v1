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
#include <stdlib.h>

struct reduced_mesh refine_reduced(
    unsigned elem_dim,
    unsigned nelem,
    unsigned nvert,
    unsigned const* elem_verts,
    double const* coords,
    double (*size_function)(double const x[]))
{
  double* sizes = malloc(sizeof(double) * nvert);
  for (unsigned i = 0; i < nvert; ++i)
    sizes[i] = size_function(coords + i * 3);
  unsigned nedges;
  unsigned* verts_of_edges;
  derive_edges(
      elem_dim,
      nelem,
      nvert,
      elem_verts,
      &nedges,
      &verts_of_edges);
  double* edge_sizes = measure_edges(
      nedges,
      verts_of_edges,
      coords,
      sizes);
  free(sizes);
  unsigned* candidates = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i)
    candidates[i] = edge_sizes[i] > 1.0;
  unsigned* edges_of_verts_offsets;
  unsigned* edges_of_verts;
  up_from_down(
      1,
      0,
      nedges,
      nvert,
      verts_of_edges,
      &edges_of_verts_offsets,
      &edges_of_verts,
      0);
  unsigned* edges_of_elems = reflect_down(
      elem_dim,
      1,
      nelem,
      nedges,
      elem_verts,
      edges_of_verts_offsets,
      edges_of_verts);
  unsigned* elems_of_edges_offsets;
  unsigned* elems_of_edges;
  up_from_down(
      elem_dim,
      1,
      nelem,
      nedges,
      edges_of_elems,
      &elems_of_edges_offsets,
      &elems_of_edges,
      0);
  unsigned* edges_of_edges_offsets;
  unsigned* edges_of_edges;
  get_star(
      1,
      elem_dim,
      nedges,
      elems_of_edges_offsets,
      elems_of_edges,
      edges_of_elems,
      &edges_of_edges_offsets,
      &edges_of_edges);
  unsigned* indset = find_indset(
      nedges,
      edges_of_edges_offsets,
      edges_of_edges,
      candidates,
      edge_sizes);
  free(edges_of_edges_offsets);
  free(edges_of_edges);
  free(candidates);
  free(edge_sizes);
  unsigned* gen_offset_of_edges = ints_exscan(indset, nedges);
  free(candidates);
  unsigned* gen_vert_of_edges = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i)
    if (gen_offset_of_edges[i] != gen_offset_of_edges[i + 1])
      gen_vert_of_edges[i] = nvert + gen_offset_of_edges[i];
  unsigned* gen_offset_of_elems;
  unsigned* gen_direction_of_elems;
  unsigned* gen_vert_of_elems;
  project_splits_to_elements(
      elem_dim,
      1,
      nelem,
      edges_of_elems,
      gen_offset_of_edges,
      gen_vert_of_edges,
      &gen_offset_of_elems,
      &gen_direction_of_elems,
      &gen_vert_of_elems);
  unsigned ngen_elems;
  unsigned* verts_of_gen_elems;
  refine_topology(
      elem_dim,
      1,
      elem_dim,
      nelem,
      elem_verts,
      gen_offset_of_elems,
      gen_vert_of_elems,
      gen_direction_of_elems,
      &ngen_elems,
      &verts_of_gen_elems);
  free(gen_offset_of_elems);
  free(gen_vert_of_elems);
  free(gen_direction_of_elems);
  double* gen_coords = refine_nodal(
      1,
      nedges,
      verts_of_edges,
      gen_offset_of_edges,
      3,
      coords);
  (void)gen_coords;
  struct reduced_mesh rm;
  return rm;
}
