#ifndef MESH_DIFF_H
#define MESH_DIFF_H

struct mesh;

unsigned mesh_diff(struct mesh* a, struct mesh* b,
    double tol, double floor);

#endif
