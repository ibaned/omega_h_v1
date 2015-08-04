#include "refine_reduced.h"
#include "derive_edges.h"
#include "measure_edges.h"
#include "refine_nodal.h"
#include "up_from_down.h"
#include "star.h"
#include "reflect_down.h"
#include "independent_set.h"
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
  struct star edge_edges = get_star(
      1,
      elem_dim,
      nedges,
      elems_of_edges_offsets,
      elems_of_edges,
      edges_of_elems);
  unsigned* indset = find_independent_set(
      nedges,
      edge_edges.offsets,
      edge_edges.edges,
      candidates,
      edge_sizes);
  free(edge_edges.offsets);
  free(edge_edges.edges);
  free(candidates);
  free(edge_sizes);
  unsigned* split_offsets = ints_exscan(indset, nedges);
  free(candidates);
  unsigned* split_new_verts = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < nedges; ++i)
    if (split_offsets[i] != split_offsets[i + 1])
      split_new_verts[i] = nvert + split_offsets[i];
  struct splits_to_elements s2e = project_splits_to_elements(
      elem_dim,
      1,
      nelem,
      edges_of_elems,
      split_offsets,
      split_new_verts);
  struct refined_topology rt = refine_topology(
      elem_dim,
      1,
      elem_dim,
      nelem,
      elem_verts,
      s2e.elem_split_offset,
      s2e.elem_split_vert,
      s2e.elem_split_direction);
  free(s2e.elem_split_offset);
  free(s2e.elem_split_vert);
  free(s2e.elem_split_direction);
  double* gen_coords = refine_nodal(
      1,
      nedges,
      verts_of_edges,
      split_offsets,
      3,
      coords);
  (void)rt;
  (void)gen_coords;
  struct reduced_mesh rm;
  return rm;
}
