#ifndef SPLIT_SLIVERS_H
#define SPLIT_SLIVERS_H

#include "mesh.h"

unsigned split_slivers(
    unsigned sliver_dim,
    struct mesh** p_m,
    double qual_floor,
    double edge_ratio_floor);

#endif
