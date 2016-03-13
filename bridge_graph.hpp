#ifndef BRIDGE_GRAPH_HPP
#define BRIDGE_GRAPH_HPP

namespace omega_h {

void bridge_graph(
    unsigned nverts,
    unsigned const* adj_offsets,
    unsigned const* adj,
    unsigned* nedges_out,
    unsigned** verts_of_edges_out);

void bridge_dual_graph(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* elems_of_elems,
    unsigned* nsides_out,
    unsigned** elems_of_sides_out,
    unsigned** elem_side_of_sides_out);

}

#endif
