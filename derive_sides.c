#include "derive_sides.h"

#include <assert.h>

#include "loop.h"
#include "tables.h"

unsigned* derive_sides(
    unsigned elem_dim,
    unsigned nsides,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_sides,
    unsigned const* elem_side_of_sides)
{
  assert(2 <= elem_dim);
  unsigned nverts_per_elem = the_down_degrees[elem_dim][0];
  unsigned nverts_per_side = the_down_degrees[elem_dim - 1][0];
  unsigned* verts_of_sides = LOOP_MALLOC(unsigned, nsides * nverts_per_side);
  unsigned const* const* elem_verts_of_sides =
    the_canonical_orders[elem_dim][elem_dim - 1][0];
  for (unsigned i = 0; i < nsides; ++i) {
    unsigned elem = elems_of_sides[i * 2];
    unsigned const* verts_of_elem = verts_of_elems + elem * nverts_per_elem;
    unsigned elem_side = elem_side_of_sides[i];
    unsigned const* elem_verts_of_side = elem_verts_of_sides[elem_side];
    unsigned* verts_of_side = verts_of_sides + i * nverts_per_side;
    for (unsigned j = 0; j < nverts_per_side; ++j) {
      unsigned elem_vert_of_side = elem_verts_of_side[j];
      unsigned vert = verts_of_elem[elem_vert_of_side];
      verts_of_side[j] = vert;
    }
  }
  return verts_of_sides;
}
