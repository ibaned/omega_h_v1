#ifndef EDGE_RING_H
#define EDGE_RING_H

#include "loop.h"

LOOP_IN unsigned find_edge_ring(
    unsigned edge,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    unsigned* edge_v,
    unsigned* ring_v);

#endif
