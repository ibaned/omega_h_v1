#ifndef REFINE_REDUCED_H
#define REFINE_REDUCED_H

#include "rv_mesh.h"

struct rv_mesh refine_reduced(
    struct rv_mesh in,
    double (*size_function)(double const x[]));

#endif
