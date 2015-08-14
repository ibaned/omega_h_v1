#ifndef BRIDGE_GRAPH_H
#define BRIDGE_GRAPH_H

void bridge_graph(
    unsigned nverts,
    unsigned const adj_offsets[],
    unsigned const adj[],
    unsigned* nedges_out,
    unsigned** verts_of_edges_out);

void bridge_dual_graph(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const elems_to_elems[],
    unsigned* nfaces_out,
    unsigned** elems_of_faces_out,
    unsigned** elems_of_faces_directions_out);

#endif
