#ifndef DERIVE_FACES_H
#define DERIVE_FACES_H

unsigned* derive_sides(
    unsigned elem_dim,
    unsigned nsides,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_sides,
    unsigned const* elem_side_of_sides);

#endif
