#ifndef REFINE_FIT_H
#define REFINE_FIT_H

struct mesh;

void refine_fit(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

#endif

