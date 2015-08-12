#include "refine_reduced.h"
#include "derive_edges.h"
#include "measure_edges.h"
#include "refine_common.h"
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
  unsigned nedges;
  unsigned* verts_of_edges;
  derive_edges(elem_dim, nelems, nverts, verts_of_elems,
      &nedges, &verts_of_edges);
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
