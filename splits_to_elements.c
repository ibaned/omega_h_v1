#include "splits_to_elements.h"
#include "tables.h"
#include "ints.h"
#include <stdlib.h>
#include <assert.h>

struct splits_to_elements project_splits_to_elements(
    unsigned elem_dim,
    unsigned split_dim,
    unsigned nelem,
    unsigned const* elem_splits,
    unsigned const* split_offsets,
    unsigned const* split_new_verts)
{
  assert(elem_dim >= split_dim);
  unsigned nelem_split = the_down_degrees[elem_dim][split_dim];
  struct splits_to_elements out;
  unsigned* elem_will_split = malloc(sizeof(unsigned) * nelem);
  ints_zero(elem_will_split, nelem);
  out.elem_split_vert = malloc(sizeof(unsigned) * nelem);
  out.elem_split_direction = malloc(sizeof(unsigned) * nelem);
  for (unsigned i = 0; i < nelem; ++i) {
    unsigned const* elem_split = elem_splits + i * nelem_split;
    for (unsigned j = 0; j < nelem_split; ++j) {
      unsigned split = elem_split[j];
      if (split_offsets[split] == split_offsets[split + 1])
        continue;
      elem_will_split[i] = 1;
      out.elem_split_vert[i] = split_new_verts[split];
      out.elem_split_direction[i] = j;
      break;
    }
  }
  out.elem_split_offset = ints_exscan(elem_will_split, nelem);
  free(elem_will_split);
  return out;
}
