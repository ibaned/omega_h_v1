#ifndef MESH_DIFF_HPP
#define MESH_DIFF_HPP

namespace omega_h {

struct mesh;

unsigned mesh_diff(struct mesh* a, struct mesh* b,
    double tol, double floor, unsigned allow_superset);

}

#endif
