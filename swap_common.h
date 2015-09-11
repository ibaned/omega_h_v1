#ifndef SWAP_COMMON_H
#define SWAP_COMMON_H

struct mesh;

unsigned swap_common(
    struct mesh** p_m,
    unsigned* candidates,
    double good_qual,
    double valid_qual,
    unsigned require_better);

#endif
