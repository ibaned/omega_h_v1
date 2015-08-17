#ifndef REFINE_BY_SIZE_H
#define REFINE_BY_SIZE_H

#include "mesh.h"

unsigned refine_by_size(
    struct mesh** p_m,
    double (*size_function)(double const x[]));

#endif
