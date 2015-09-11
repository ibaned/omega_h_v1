#ifndef SWAP_TOPOLOGY_H
#define SWAP_TOPOLOGY_H

unsigned* swap_topology(
    unsigned nedges,
    unsigned const* candidates,
    unsigned const* gen_offset_of_edges,
    unsigned const* edge_codes,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets);

#endif
