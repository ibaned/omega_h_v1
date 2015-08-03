#include "dual.h"
#include "tables.h"
#include "intersect.h"
#include <assert.h>
#include <stdlib.h>

/* this is very similar to reflect_down */

unsigned* get_dual_using_verts(
    unsigned elem_dim,
    unsigned nelem,
    unsigned const* elem_verts,
    unsigned const* vert_elem_offsets,
    unsigned const* vert_elems)
{
  assert(elem_dim > 0);
  unsigned side_dim = elem_dim - 1;
  unsigned nelem_vert = the_down_degrees[elem_dim][0];
  unsigned nelem_side = the_down_degrees[elem_dim][side_dim];
  unsigned nside_vert = the_down_degrees[side_dim][0];
  unsigned const* const* any_side_vert =
    the_canonical_orders[elem_dim][side_dim][0];
  unsigned* out = malloc(sizeof(unsigned) * nelem * nelem_side);
  for (unsigned i = 0; i < nelem; ++i) {
    unsigned const* elem_vert = elem_verts + i * nelem_vert;
    unsigned* elem_out = out + i * nelem_side;
    for (unsigned j = 0; j < nelem_side; ++j) {
      unsigned const* side_vert = any_side_vert[j];
      unsigned elem_buf[MAX_UP];
      unsigned nelem_buf = 0;
      for (unsigned k = 0; k < nside_vert; ++k) {
        unsigned vert = elem_vert[side_vert[k]];
        unsigned first_elem = vert_elem_offsets[vert];
        unsigned end_elem = vert_elem_offsets[vert + 1];
        if (nelem_buf)
          nelem_buf = intersect(elem_buf, nelem_buf,
              vert_elems + first_elem,
              end_elem - first_elem);
        else
          nelem_buf = copy(
              vert_elems + first_elem,
              elem_buf,
              end_elem - first_elem);
        assert(nelem_buf);
      }
      assert(nelem_buf <= 2);
      if (nelem_buf == 1)
        elem_out[i] = DUAL_NONE;
      else
        elem_out[i] = (elem_buf[0] == i) ? elem_buf[1] : elem_buf[0];
    }
  }
  return out;
}
