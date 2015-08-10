#include "verify.h"
#include "element_qualities.h"
#include "tables.h"
#include "doubles.h"
#include <stdlib.h>
#include <assert.h>

#define EPSILON 1e-10

void verify(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords)
{
  assert(elem_dim >= 2);
  assert(elem_dim <= 3);
  assert(nelems);
  assert(nverts);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  for (unsigned i = 0; i < nelems * verts_per_elem; ++i)
    assert(verts_of_elems[i] < nverts);
  for (unsigned i = 0; i < nverts * 3; ++i) {
    assert(coords[i] < (1 + EPSILON));
    assert(coords[i] > -EPSILON);
  }
  double minq = min_element_quality(elem_dim, nelems, verts_of_elems, coords);
  assert(minq > 0);
}
