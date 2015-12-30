#include "collapses_to_elements.h"

#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

void collapses_to_elements(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned** p_gen_offset_of_elems,
    unsigned** p_gen_vert_of_elems,
    unsigned** p_gen_direction_of_elems,
    unsigned** p_offset_of_same_elems,
    unsigned** p_dead_sides)
{
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned const* sides_of_elems = 0;
  unsigned* dead_sides = 0;
  if (p_dead_sides) {
    sides_of_elems = mesh_ask_down(m, elem_dim, elem_dim - 1);
    unsigned nsides = mesh_count(m, elem_dim - 1);
    dead_sides = uints_filled(nsides, 0);
  }
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned sides_per_elem = the_down_degrees[elem_dim][elem_dim - 1];
  unsigned* elem_will_gen = uints_filled(nelems, 0);
  unsigned* gen_vert_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* gen_direction_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* elem_is_same = LOOP_MALLOC(unsigned, nelems);
  unsigned const* elem_sides_opp_verts = the_opposite_orders[elem_dim][0];
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
    unsigned gen_dir = INVALID;
    for (unsigned j = 0; j < verts_per_elem; ++j)
      if (verts_of_elem[j] == gen_vert)
        gen_dir = j;
    unsigned elem_will_die = (gen_dir != INVALID);
    if (elem_will_die) {
      if (p_dead_sides) {
        unsigned elem_side = elem_sides_opp_verts[gen_dir];
        unsigned side = sides_of_elems[i * sides_per_elem + elem_side];
        dead_sides[side] = 1;
      }
      continue;
    }
    elem_will_gen[i] = 1;
    gen_vert_of_elems[i] = gen_vert;
    gen_direction_of_elems[i] = direction;
  }
  *p_gen_offset_of_elems = uints_exscan(elem_will_gen, nelems);
  loop_free(elem_will_gen);
  *p_gen_vert_of_elems = gen_vert_of_elems;
  *p_gen_direction_of_elems = gen_direction_of_elems;
  *p_offset_of_same_elems = uints_exscan(elem_is_same, nelems);
  loop_free(elem_is_same);
  if (p_dead_sides)
    *p_dead_sides = dead_sides;
}
