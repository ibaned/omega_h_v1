#ifndef WARP_TO_LIMIT
#define WARP_TO_LIMIT

namespace omega_h {

struct mesh;

unsigned mesh_warp_to_limit(struct mesh* m, double qual_floor);

}

#endif
