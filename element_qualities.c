#include "element_qualities.h"
#include "tables.h"
#include "algebra.h"
#include "quality.h"
#include "doubles.h"
#include <stdlib.h>

double* element_qualities(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords)
{
  double* out = malloc(sizeof(double) * nelems);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  quality_function qf = the_quality_functions[elem_dim];
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    double elem_x[MAX_DOWN][3];
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned vert = verts_of_elem[j];
      copy_vector(coords + vert * 3, elem_x[j], 3);
    }
    out[i] = qf(elem_x);
  }
  return out;
}

double min_element_quality(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords)
{
  double* quals = element_qualities(elem_dim, nelems, verts_of_elems, coords);
  double mq = doubles_min(quals, nelems);
  free(quals);
  return mq;
}
