#ifndef MARK_H
#define MARK_H

unsigned* mark_down(
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* marked_highs);

unsigned* mark_up(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* lows_of_highs,
    unsigned const* marked_lows);

struct mesh;

unsigned* mesh_mark_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* marked_highs);
unsigned* mesh_mark_up(struct mesh* m, unsigned low_dim, unsigned high_dim,
    unsigned const* marked_lows);

void mesh_mark_dual_layers(
    struct mesh* m,
    unsigned** marked,
    unsigned nlayers);

void unmark_boundary(
    unsigned elem_dim,
    unsigned ent_dim,
    unsigned nents,
    unsigned const* verts_of_ents,
    unsigned const* vert_class_dim,
    unsigned* marked);

#endif
