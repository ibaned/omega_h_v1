#ifndef REFINE_CLASS_H
#define REFINE_CLASS_H

struct mesh;

void refine_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned* offset_of_doms[4],
    unsigned ngen_ents[4][4],
    unsigned gen_offsets[4][5]);

#endif
