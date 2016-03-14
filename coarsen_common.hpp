#ifndef COARSEN_COMMON_HPP
#define COARSEN_COMMON_HPP

namespace omega_h {

struct mesh;

unsigned coarsen_common(
    struct mesh* m,
    double quality_floor,
    unsigned require_better);

}

#endif
