#ifndef REFINE_COMMON_HPP
#define REFINE_COMMON_HPP

struct mesh;

unsigned refine_common(
    struct mesh* m,
    unsigned src_dim,
    double qual_floor,
    unsigned require_better);

#endif
