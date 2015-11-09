#ifndef OWNERS_FROM_VERTS_H
#define OWNERS_FROM_VERTS_H

void owners_from_verts(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    unsigned const* own_rank_of_verts,
    unsigned const* own_idx_of_verts,
    unsigned** p_own_ranks,
    unsigned** p_own_idxs);

#endif
