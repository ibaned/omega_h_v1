#ifndef WARP_TO_LIMIT
#define WARP_TO_LIMIT

double* warp_to_limit(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords,
    double const* warps);

#endif
