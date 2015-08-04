#ifndef REFINE_NODAL_H
#define REFINE_NODAL_H

double* refine_nodal(
    unsigned split_dim,
    unsigned nsplit,
    unsigned const* split_verts,
    unsigned const* split_offsets,
    unsigned ncomp,
    double const* field);

#endif
