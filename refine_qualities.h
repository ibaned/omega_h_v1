#ifndef REFINE_QUALITIES_H
#define REFINE_QUALITIES_H

double* refine_qualities(
    unsigned elem_dim,
    unsigned split_dim,
    unsigned nsplit,
    unsigned const* split_verts,
    unsigned const* split_elem_offsets,
    unsigned const* split_elem_edges,
    unsigned const* split_elem_directions,
    unsigned const* elem_verts,
    double const* coords);

#endif
