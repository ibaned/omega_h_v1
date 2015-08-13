#ifndef SPLIT_SLIVER_TRIS
#define SPLIT_SLIVER_TRIS

int split_sliver_tris(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    unsigned** p_class_dim,
    double qual_floor,
    double edge_ratio_floor);

#endif
