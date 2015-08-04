#ifndef BRIDGE_GRAPH_H
#define BRIDGE_GRAPH_H

void bridge_graph(
    unsigned nverts,
    unsigned const adj_offsets[],
    unsigned const adj[],
    unsigned* nedges_out,
    unsigned** verts_of_edges_out);

#endif
