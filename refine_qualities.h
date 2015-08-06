#ifndef REFINE_QUALITIES_H
#define REFINE_QUALITIES_H

double* refine_qualities(
    unsigned elem_dim,
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_srcs_offsets,
    unsigned const* elems_of_srcs,
    unsigned const* elems_of_srcs_directions,
    unsigned const* candidate_srcs,
    double const* coords);

#endif
