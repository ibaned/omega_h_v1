#ifndef OWNERS_FROM_VERTS_H
#define OWNERS_FROM_VERTS_H

void owners_from_verts(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* own_part_of_verts,
    unsigned const* own_idx_of_verts,
    unsigned** p_own_parts,
    unsigned** p_own_idxs);

#endif
