#ifndef SHUFFLE_MESH_H
#define SHUFFLE_MESH_H

struct mesh;

void shuffle_mesh(struct mesh* m, unsigned const* old_to_new_verts);

#endif
