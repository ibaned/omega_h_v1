#ifndef SWAP_TOPOLOGY_HPP
#define SWAP_TOPOLOGY_HPP

unsigned* get_swap_topology_offsets(
    unsigned ent_dim,
    unsigned nedges,
    unsigned const* indset,
    unsigned const* ring_sizes);

struct mesh;

unsigned* mesh_swap_topology(
    struct mesh* m,
    unsigned ent_dim,
    unsigned const* candidates,
    unsigned const* gen_offset_of_edges);

#endif
