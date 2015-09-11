#ifndef SWAP_QUALITIES_H
#define SWAP_QUALITIES_H

void swap_qualities(
    unsigned nedges,
    unsigned* candidates,
    unsigned const* tets_of_edges_offsets,
    unsigned const* tets_of_edges,
    unsigned const* tets_of_edges_directions,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_tets,
    double const* coords,
    double good_qual,
    double valid_qual,
    double const* elem_quals,
    unsigned require_better,
    double** p_qualities,
    unsigned** p_codes,
    unsigned** p_gen_elems_per_edge);

#endif
