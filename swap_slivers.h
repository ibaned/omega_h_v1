#ifndef SWAP_SLIVERS_H
#define SWAP_SLIVERS_H

struct mesh;

unsigned swap_slivers(
    struct mesh** p_m,
    double good_qual,
    unsigned nlayers);

#endif
