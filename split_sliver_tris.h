#ifndef SPLIT_SLIVER_TRIS
#define SPLIT_SLIVER_TRIS

#include "mesh.h"

unsigned split_sliver_tris(
    struct mesh** p_m,
    double qual_floor,
    double edge_ratio_floor);

#endif
