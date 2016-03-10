#ifndef MARK_HPP
#define MARK_HPP

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

unsigned* mesh_mark_down_local(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* marked_highs);
unsigned* mesh_mark_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* marked_highs);
unsigned* mesh_mark_up(struct mesh* m, unsigned low_dim, unsigned high_dim,
    unsigned const* marked_lows);

void mesh_mark_dual_layers(
    struct mesh* m,
    unsigned** marked,
    unsigned nlayers);

unsigned* mark_class(
    unsigned nents,
    unsigned target_dim,
    unsigned target_id,
    unsigned const* class_dim_of_ents,
    unsigned const* class_id_of_ents);

unsigned* mesh_mark_class(struct mesh* m, unsigned ent_dim,
    unsigned target_dim, unsigned target_id);

unsigned* mesh_mark_class_closure_verts(struct mesh* m, unsigned target_dim,
    unsigned target_id);

unsigned* mesh_mark_slivers(struct mesh* m, double good_qual, unsigned nlayers);

unsigned* mark_part_boundary(
    unsigned nsides,
    unsigned const* elems_of_sides_offsets);

void mesh_unmark_boundary(
    struct mesh* m,
    unsigned ent_dim,
    unsigned* marked);

unsigned* mesh_mark_part_boundary(struct mesh* m);

#endif
