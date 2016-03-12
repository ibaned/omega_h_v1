#ifndef REORDER_MESH_HPP
#define REORDER_MESH_HPP

struct mesh;

void reorder_mesh(struct mesh* m, unsigned const* old_to_new_verts);

#endif
