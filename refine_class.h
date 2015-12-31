#ifndef REFINE_CLASS_H
#define REFINE_CLASS_H

struct mesh;

void refine_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

#endif
