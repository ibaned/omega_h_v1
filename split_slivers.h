#ifndef SPLIT_SLIVERS_H
#define SPLIT_SLIVERS_H

#include "mesh.h"
#include "quality.h"

unsigned split_slivers(
    struct mesh** p_m,
    unsigned sliver_dim,
    enum sliver_type st,
    double good_qual,
    double valid_qual);

#endif
