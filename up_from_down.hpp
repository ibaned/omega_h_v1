#ifndef UP_FROM_DOWN_HPP
#define UP_FROM_DOWN_HPP

/* given a downward adjacency graph,
 * computes its inverse upward adjacency graph,
 * including information about the relative
 * orientation of entities (directions)
 */

void up_from_down(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned nlows,
    unsigned const* lows_of_highs,
    unsigned** offsets_out,
    unsigned** highs_of_lows_out,
    unsigned** directions_out);

#endif
