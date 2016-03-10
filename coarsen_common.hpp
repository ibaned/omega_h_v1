#ifndef COARSEN_COMMON_HPP
#define COARSEN_COMMON_HPP

struct mesh;

unsigned coarsen_common(
    struct mesh* m,
    double quality_floor,
    unsigned require_better);

#endif
