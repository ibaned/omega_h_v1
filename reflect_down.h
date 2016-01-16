#ifndef REFLECT_DOWN_H
#define REFLECT_DOWN_H

/* given two entity types (high) and (low), (high >= low)
 * as well as the adjacencies from (high) to vertices
 * and from vertices to (low),
 * this function derives the downward adjacency from (high) to (low).
 */

unsigned* reflect_down(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* verts_of_highs,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts);

struct mesh;

unsigned* mesh_reflect_down(
    struct mesh* m,
    unsigned high_dim,
    unsigned low_dim);

#endif
