#ifndef SMOOTH_HPP
#define SMOOTH_HPP

namespace omega_h {

struct mesh;

unsigned mesh_smooth_field(struct mesh* m, char const* name,
    double tol, unsigned maxiter);

}

#endif
