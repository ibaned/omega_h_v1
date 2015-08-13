#ifndef REFINE_COMMON_H
#define REFINE_COMMON_H

void refine_common(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    unsigned** p_class_dim,
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* candidates,
    unsigned const* srcs_of_elems,
    unsigned const* elems_of_srcs_offsets,
    unsigned const* elems_of_srcs,
    unsigned const* elems_of_srcs_directions,
    unsigned const* srcs_of_srcs_offsets,
    unsigned const* srcs_of_srcs);

#endif
