#ifndef COARSEN_CONSERVE_HPP
#define COARSEN_CONSERVE_HPP

struct mesh;

void coarsen_conserve(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_offset_of_elems,
    unsigned const* offset_of_same_elems);

#endif
