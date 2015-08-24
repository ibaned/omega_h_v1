#ifndef WARP_TO_LIMIT
#define WARP_TO_LIMIT

void warp_to_limit(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords,
    double const* warps,
    double** p_coords,
    double** p_warps);

struct mesh;

void mesh_warp_to_limit(struct mesh* m);

#endif
