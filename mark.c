#include "mark.h"
#include "ints.h"
#include "mesh.h"
#include <stdlib.h>

unsigned* mark_down(
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* marked_highs)
{
  unsigned* low_marks = malloc(sizeof(unsigned) * nlows);
  for (unsigned i = 0; i < nlows; ++i) {
    unsigned first_use = highs_of_lows_offsets[i];
    unsigned end_use = highs_of_lows_offsets[i + 1];
    low_marks[i] = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned high = highs_of_lows[j];
      if (marked_highs[high]) {
        low_marks[i] = 1;
        break;
      }
    }
  }
  return low_marks;
}

unsigned* mesh_mark_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* marked_highs)
{
  return mark_down(mesh_count(m, low_dim),
      mesh_ask_up(m, low_dim, high_dim)->offsets,
      mesh_ask_up(m, low_dim, high_dim)->adj,
      marked_highs);
}
