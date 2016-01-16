#ifndef DUAL_H
#define DUAL_H

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

struct mesh;

unsigned* mesh_get_dual(struct mesh* m);

#endif
