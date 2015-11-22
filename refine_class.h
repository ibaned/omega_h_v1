#ifndef REFINE_CLASS_H
#define REFINE_CLASS_H

struct mesh;

void refine_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned const* gen_offset_of_srcs);

#endif
