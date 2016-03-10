#ifndef STAR_HPP
#define STAR_HPP

/* given a dual pair of upward and downward
 * adjacencies between two dimensions,
 * this function returns the "second order adjacency"
 * of "low" entities adjacent to one another
 * via a "high" entity.
 *
 * the term star comes from applying this function
 * using mesh vertices and edges, in which case one gets
 * at each vertex a "star" of vertices adjacent
 * across an edge.
 */

struct mesh;

void mesh_get_star(
    struct mesh* m,
    unsigned low_dim,
    unsigned high_dim,
    unsigned** p_star_offsets,
    unsigned** p_star);

#endif
