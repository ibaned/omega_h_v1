#ifndef EDGE_RING_HPP
#define EDGE_RING_HPP

#include "loop.hpp"

namespace omega_h {

LOOP_IN unsigned find_edge_ring(
    unsigned edge,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    unsigned* edge_v,
    unsigned* ring_v);

}

#endif
