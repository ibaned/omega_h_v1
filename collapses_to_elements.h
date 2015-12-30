#ifndef COLLAPSES_TO_ELEMENTS_H
#define COLLAPSES_TO_ELEMENTS_H

struct mesh;

void collapses_to_elements(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned** gen_offset_of_elems_out,
    unsigned** gen_vert_of_elems_out,
    unsigned** gen_direction_of_elems_out,
    unsigned** offset_of_same_elems_out,
    unsigned** p_dead_sides);

#endif
