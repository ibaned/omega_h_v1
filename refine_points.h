#ifndef REFINE_POINTS_H
#define REFINE_POINTS_H

void refine_points(
    unsigned elem_dim,
    unsigned nsrc_elems,
    unsigned nnew_elems,
    unsigned const* verts_of_new_elems,
    double const* coords_of_verts,
    unsigned const* pts_of_src_elems_offsets,
    unsigned const* pts_of_src_elems,
    double const* coords_of_pts,
    unsigned** p_pts_of_new_elems_offsets,
    unsigned** p_pts_of_new_elems);

#endif
