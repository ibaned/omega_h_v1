#include "collapses_to_elements.h"
#include "tables.h"
#include "ints.h"
#include <stdlib.h>

void collapses_to_elements(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned** gen_offset_of_elems_out,
    unsigned** gen_vert_of_elems_out,
    unsigned** gen_direction_of_elems_out,
    unsigned** offset_of_same_elems_out)
{
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* elem_will_gen = malloc(sizeof(unsigned) * nelems);
  unsigned* gen_vert_of_elems = malloc(sizeof(unsigned) * nelems);
  unsigned* gen_direction_of_elems = malloc(sizeof(unsigned) * nelems);
  unsigned* elem_is_same = malloc(sizeof(unsigned) * nelems);
  ints_zero(elem_will_gen, nelems);
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    unsigned col_vert = INVALID;
    unsigned direction = INVALID;
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned vert = verts_of_elem[j];
      if (gen_offset_of_verts[vert] != gen_offset_of_verts[vert + 1]) {
        col_vert = vert;
        direction = j;
      }
    }
    if (direction == INVALID) {
      elem_is_same[i] = 1;
      /* not adjacent to collapsing vertex */
      continue;
    }
    elem_is_same[i] = 0;
    unsigned gen_vert = gen_vert_of_verts[col_vert];
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      if (verts_of_elem[j] == gen_vert)
        /* has both vertices of the collapsing edge,
           will be destroyed by the collapse */
        continue;
    }
    elem_will_gen[i] = 1;
    gen_vert_of_elems[i] = gen_vert;
    gen_direction_of_elems[i] = direction;
  }
  *gen_offset_of_elems_out = ints_exscan(elem_will_gen, nelems);
  free(elem_will_gen);
  *gen_vert_of_elems_out = gen_vert_of_elems;
  *gen_direction_of_elems_out = gen_direction_of_elems;
  *offset_of_same_elems_out = ints_exscan(elem_is_same, nelems);
  free(elem_is_same);
}
