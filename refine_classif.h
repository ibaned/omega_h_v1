#ifndef REFINE_CLASSIF_H
#define REFINE_CLASSIF_H

unsigned* refine_classif(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned const* class_dim);

#endif
