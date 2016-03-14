#ifndef SWAP_HPP
#define SWAP_HPP

namespace omega_h {

struct mesh;

unsigned swap_slivers(
    struct mesh* m,
    double good_qual,
    unsigned nlayers);

}

#endif
