#ifndef SWAP_CONSERVE_HPP
#define SWAP_CONSERVE_HPP

namespace omega_h {

struct mesh;

void swap_conserve(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_edges,
    unsigned const* offset_of_same_elems);

}

#endif
