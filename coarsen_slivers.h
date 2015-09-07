#ifndef COARSEN_SLIVERS_H
#define COARSEN_SLIVERS_H

struct mesh;

unsigned coarsen_slivers(
    struct mesh** p_m,
    double quality_floor);

#endif
