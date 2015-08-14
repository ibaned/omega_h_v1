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

/* compute the dual or element-to-element graph.
 * this may not find an element in a given direction,
 * so those entries will be marked INVALID.
 */

unsigned* get_dual(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts);

#endif
