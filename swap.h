#ifndef SWAP_H
#define SWAP_H

struct mesh;

unsigned swap_common(
    struct mesh** p_m,
    unsigned* candidates);

unsigned swap_slivers(
    struct mesh** p_m,
    double good_qual,
    unsigned nlayers);

#endif
