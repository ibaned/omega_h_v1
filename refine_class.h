#ifndef REFINE_CLASS_H
#define REFINE_CLASS_H

unsigned* refine_class(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned const* class_dim);

#endif
