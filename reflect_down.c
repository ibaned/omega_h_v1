#include "reflect_down.h"
#include "intersect.h"
#include "tables.h"
#include <stdlib.h>
#include <assert.h>

unsigned* reflect_down(
    unsigned up_dim,
    unsigned down_dim,
    unsigned nup,
    unsigned ndown,
    unsigned const* up_verts,
    unsigned const* vert_down_offsets,
    unsigned const* vert_downs)
{
  unsigned nup_vert = the_down_degrees[up_dim][0];
  unsigned nup_down = the_down_degrees[up_dim][down_dim];
  unsigned ndown_vert = the_down_degrees[down_dim][0];
  unsigned* out = malloc(sizeof(unsigned) * ndown * ndown_vert);
  unsigned const* const* any_down_vert = the_canonical_orders[up_dim][down_dim][0];
  for (unsigned i = 0; i < nup; ++i) {
    unsigned const* up_vert = up_verts + i * nup_vert;
    unsigned* up_out = out + i * nup_down;
    for (unsigned j = 0; j < nup_down; ++j) {
      unsigned const* down_vert = any_down_vert[j];
      unsigned up_buf[MAX_UP];
      unsigned nup_buf = 0;
      for (unsigned k = 0; k < ndown_vert; ++k) {
        unsigned vert = up_vert[down_vert[k]];
        unsigned first_down = vert_down_offsets[vert];
        unsigned end_down = vert_down_offsets[vert + 1];
        if (nup_buf)
          nup_buf = intersect(up_buf, nup_buf,
              vert_downs + first_down,
              end_down - first_down);
        else
          nup_buf = copy(
              vert_downs + first_down,
              up_buf,
              end_down - first_down);
        assert(nup_buf);
      }
      assert(nup_buf == 1);
      up_out[j] = up_buf[0];
    }
  }
  return out;
}
