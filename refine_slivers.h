#ifndef REFINE_SLIVERS_H
#define REFINE_SLIVERS_H

struct mesh;

unsigned refine_slivers(
    struct mesh** p_m,
    double good_qual,
    double valid_qual);

#endif
