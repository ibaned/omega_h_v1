#ifndef SMOOTH_HPP
#define SMOOTH_HPP

struct mesh;

unsigned mesh_smooth_field(struct mesh* m, char const* name,
    double tol, unsigned maxiter);

#endif
