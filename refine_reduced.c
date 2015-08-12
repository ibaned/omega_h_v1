#include "refine_reduced.h"
#include "refine_common.h"
#include "up_from_down.h"
#include "measure_edges.h"
#include "bridge_graph.h"
#include "star.h"
#include "ints.h"
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
  refine_common(elem_dim, p_nelems, p_nverts, p_verts_of_elems, p_coords,
      1, nedges, verts_of_edges, candidates);
  free(verts_of_edges);
  free(candidates);
  return 1;
}
