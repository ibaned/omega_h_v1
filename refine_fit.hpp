#ifndef REFINE_FIT_HPP
#define REFINE_FIT_HPP

struct mesh;

void refine_fit(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

#endif

