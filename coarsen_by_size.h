#ifndef COARSEN_BY_SIZE_H
#define COARSEN_BY_SIZE_H

#include "mesh.h"

unsigned coarsen_by_size(
    struct mesh** p_m,
    double (*size_function)(double const x[]),
    double quality_floor);

#endif
