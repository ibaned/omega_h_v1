#ifndef SWAP_QUALITIES_H
#define SWAP_QUALITIES_H

struct mesh;

void mesh_swap_qualities(
    struct mesh* m,
    unsigned* candidates,
    double** p_qualities,
    unsigned** p_codes,
    unsigned** p_ring_sizes);

#endif
