#include "collapses_to_elements.h"

#include "ints.h"    // for ints_exscan, ints_zero
#include "loop.h"  // for malloc, free
#include "tables.h"  // for INVALID, the_down_degrees

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
  unsigned* elem_will_gen = loop_malloc(sizeof(unsigned) * nelems);
  unsigned* gen_vert_of_elems = loop_malloc(sizeof(unsigned) * nelems);
  unsigned* gen_direction_of_elems = loop_malloc(sizeof(unsigned) * nelems);
  unsigned* elem_is_same = loop_malloc(sizeof(unsigned) * nelems);
  uints_zero(elem_will_gen, nelems);
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
    unsigned elem_will_die = 0;
    for (unsigned j = 0; j < verts_per_elem; ++j)
      if (verts_of_elem[j] == gen_vert)
        elem_will_die = 1;
    if (elem_will_die)
      continue;
    elem_will_gen[i] = 1;
    gen_vert_of_elems[i] = gen_vert;
    gen_direction_of_elems[i] = direction;
  }
  *gen_offset_of_elems_out = uints_exscan(elem_will_gen, nelems);
  loop_free(elem_will_gen);
  *gen_vert_of_elems_out = gen_vert_of_elems;
  *gen_direction_of_elems_out = gen_direction_of_elems;
  *offset_of_same_elems_out = uints_exscan(elem_is_same, nelems);
  loop_free(elem_is_same);
}
