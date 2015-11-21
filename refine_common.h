#ifndef REFINE_COMMON_H
#define REFINE_COMMON_H

struct mesh;

struct mesh* refine_common(
    struct mesh* m,
    unsigned src_dim,
    unsigned* candidates,
    double qual_floor,
    unsigned require_better);

#endif
