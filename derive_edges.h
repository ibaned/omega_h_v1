#ifndef DERIVE_EDGES_H
#define DERIVE_EDGES_H

void derive_edges(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    unsigned* nedges_out,
    unsigned** verts_of_edges_out);

#endif
