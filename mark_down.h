#ifndef MARK_DOWN_H
#define MARK_DOWN_H

unsigned* mark_down(
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* high_offsets);

struct mesh;

unsigned* mesh_mark_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* high_offsets);

#endif
