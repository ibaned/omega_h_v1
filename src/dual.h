#ifndef DUAL_H
#define DUAL_H

#define DUAL_NONE (~((unsigned)0))

unsigned* get_dual_using_verts(
    unsigned elem_dim,
    unsigned nelem,
    unsigned const* down_edges,
    unsigned const* up_offsets,
    unsigned const* up_edges);

#endif
