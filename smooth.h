#ifndef SMOOTH_H
#define SMOOTH_H

struct mesh;

unsigned mesh_smooth_field(struct mesh* m, char const* name,
    double tol, unsigned maxiter);

#endif
