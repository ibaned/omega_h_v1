#ifndef DUAL_HPP
#define DUAL_HPP

/* compute the dual or element-to-element graph.
 * this may not find an element in a given direction,
 * so those entries will be marked INVALID.
 */

namespace omega_h {

unsigned* dual_from_verts(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts);

unsigned* dual_from_sides(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* sides_of_elems,
    unsigned const* elems_of_sides_offsets,
    unsigned const* elems_of_sides);

struct mesh;

unsigned* mesh_get_dual_from_verts(struct mesh* m);

unsigned* mesh_get_dual_from_sides(struct mesh* m);

}

#endif
