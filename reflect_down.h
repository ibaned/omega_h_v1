#ifndef REFLECT_DOWN_H
#define REFLECT_DOWN_H

unsigned* reflect_down(
    unsigned up_dim,
    unsigned down_dim,
    unsigned nup,
    unsigned ndown,
    unsigned const* up_verts,
    unsigned const* vert_down_offsets,
    unsigned const* vert_downs);

#endif
