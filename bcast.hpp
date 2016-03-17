#ifndef BCAST_HPP
#define BCAST_HPP

namespace omega_h {

struct mesh;

struct mesh* bcast_mesh_metadata(struct mesh* m, unsigned is_source);

}

#endif
