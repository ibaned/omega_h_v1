#ifndef REFINE_CONSERVE_HPP
#define REFINE_CONSERVE_HPP

struct mesh;

void refine_conserve(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

#endif
