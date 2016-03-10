#ifndef COARSEN_FIT_H
#define COARSEN_FIT_H

struct mesh;

void coarsen_fit(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_offset_of_elems,
    unsigned const* offset_of_same_elems);

#endif

