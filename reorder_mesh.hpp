#ifndef REORDER_MESH_HPP
#define REORDER_MESH_HPP

namespace omega_h {

struct mesh;

void reorder_mesh(struct mesh* m, unsigned const* old_to_new_verts);

}

#endif
