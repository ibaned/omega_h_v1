#ifndef WARP_TO_LIMIT
#define WARP_TO_LIMIT

unsigned warp_to_limit(
    unsigned elem_dim,
    unsigned nelems,
    unsigned nverts,
    unsigned const* verts_of_elems,
    double const* coords,
    double const* warps,
    double qual_floor,
    double** p_coords,
    double** p_warps);

struct mesh;

unsigned mesh_warp_to_limit(struct mesh* m, double qual_floor);

#endif
