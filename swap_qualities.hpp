#ifndef SWAP_QUALITIES_HPP
#define SWAP_QUALITIES_HPP

namespace omega_h {

struct mesh;

void mesh_swap_qualities(
    struct mesh* m,
    unsigned* candidates,
    double** p_qualities,
    unsigned** p_ring_sizes);

}

#endif
