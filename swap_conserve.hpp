#ifndef SWAP_CONSERVE_H
#define SWAP_CONSERVE_H

struct mesh;

void swap_conserve(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_edges,
    unsigned const* offset_of_same_elems);

#endif
