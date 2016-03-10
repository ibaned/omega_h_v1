#ifndef SWAP_FIT_H
#define SWAP_FIT_H

struct mesh;

void swap_fit(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_edges,
    unsigned const* offset_of_same_elems);

#endif

