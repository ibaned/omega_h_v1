#include "mark_down.h"
#include "ints.h"
#include "mesh.h"
#include <stdlib.h>

unsigned* mark_down(
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* high_offsets)
{
  unsigned* low_marks = malloc(sizeof(unsigned) * nlows);
  for (unsigned i = 0; i < nlows; ++i) {
    unsigned first_use = highs_of_lows_offsets[i];
    unsigned end_use = highs_of_lows_offsets[i + 1];
    low_marks[i] = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned high = highs_of_lows[j];
      if (high_offsets[high] != high_offsets[high + 1]) {
        low_marks[i] = 1;
        break;
      }
    }
  }
  unsigned* out = ints_exscan(low_marks, nlows);
  free(low_marks);
  return out;
}

unsigned* mesh_mark_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* high_offsets)
{
  return mark_down(mesh_count(m, low_dim),
      mesh_ask_up(m, low_dim, high_dim)->offsets,
      mesh_ask_up(m, low_dim, high_dim)->adj,
      high_offsets);
}
