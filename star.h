#ifndef STAR_H
#define STAR_H

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

void get_star(
    unsigned low_dim,
    unsigned high_dim,
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* lows_of_highs,
    unsigned** star_offsets_out,
    unsigned** star_out);

#endif
