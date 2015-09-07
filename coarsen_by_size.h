#ifndef COARSEN_BY_SIZE_H
#define COARSEN_BY_SIZE_H

struct mesh;

unsigned coarsen_by_size(
    struct mesh** p_m,
    double quality_floor,
    double size_ratio_floor,
    unsigned require_better);

#endif
