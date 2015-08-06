#ifndef REFLECT_DOWN_H
#define REFLECT_DOWN_H

/* given two entity types (high) and (low), (high >= low)
 * as well as the adjacencies from (high) to vertices
 * and from vertices to (low),
 * this function derives the downward adjacency from (high) to (low).
 *
 * it will also handle the special case (high == low == element),
 * in which case it is computing the dual or element-to-element graph,
 * and it may not find an element in a given direction,
 * so those entries will be marked INVALID.
 */

unsigned* reflect_down(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned nlows,
    unsigned const* verts_of_highs,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts);

#endif
