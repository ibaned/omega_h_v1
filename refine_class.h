#ifndef REFINE_CLASS_H
#define REFINE_CLASS_H

void refine_class(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned const* dims_in,
    unsigned const* ids_in,
    unsigned** p_dims_out,
    unsigned** p_ids_out);

#endif
