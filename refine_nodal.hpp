#ifndef REFINE_NODAL_H
#define REFINE_NODAL_H

double* refine_nodal(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned comps_per_vert,
    double const* field);

#endif
