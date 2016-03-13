#ifndef COARSEN_HPP
#define COARSEN_HPP

namespace omega_h {

struct mesh;

unsigned coarsen_by_size(
    struct mesh* m,
    double quality_floor,
    double size_ratio_floor);

unsigned coarsen_slivers(
    struct mesh* m,
    double quality_floor,
    unsigned nlayers);

}

#endif
