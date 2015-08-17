#ifndef REFINE_COMMON_H
#define REFINE_COMMON_H

#include "mesh.h"

void refine_common(
    struct mesh** p_m,
    unsigned src_dim,
    unsigned const* candidates);

#endif
