#ifndef REFLECT_DOWN_HPP
#define REFLECT_DOWN_HPP

namespace omega_h {

/* given two entity types (high) and (low), (high >= low)
 * as well as the adjacencies from (high) to vertices
 * and from vertices to (low),
 * this function derives the downward adjacency from (high) to (low).
 */

struct mesh;

unsigned* mesh_reflect_down(
    struct mesh* m,
    unsigned high_dim,
    unsigned low_dim);

}

#endif
