#ifndef REFINE_COMMON_H
#define REFINE_COMMON_H

struct mesh;

unsigned refine_common(
    struct mesh** p_m,
    unsigned src_dim,
    unsigned const* candidates,
    double qual_floor,
    unsigned require_better);

#endif
