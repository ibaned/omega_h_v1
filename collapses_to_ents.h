#ifndef COLLAPSES_TO_ENTS_H
#define COLLAPSES_TO_ENTS_H

struct mesh;

void collapses_to_ents(
    struct mesh* m,
    unsigned ent_dim,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned const* dead_ents,
    unsigned** gen_offset_of_ents_out,
    unsigned** gen_vert_of_ents_out,
    unsigned** gen_direction_of_ents_out,
    unsigned** offset_of_same_ents_out,
    unsigned** p_dead_sides);

#endif
