#ifndef REFINE_SLIVERS_H
#define REFINE_SLIVERS_H

struct mesh;

unsigned refine_slivers(
    struct mesh** p_m,
    unsigned src_dim,
    double good_qual,
    double valid_qual,
    unsigned require_better);

#endif
