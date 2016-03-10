#ifndef REFINE_COMMON_H
#define REFINE_COMMON_H

struct mesh;

unsigned refine_common(
    struct mesh* m,
    unsigned src_dim,
    double qual_floor,
    unsigned require_better);

#endif
