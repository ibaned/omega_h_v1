#ifndef GHOST_MESH_H
#define GHOST_MESH_H

struct mesh;

void ghost_mesh(struct mesh* m, unsigned nlayers);
void unghost_mesh(struct mesh* m);

void mesh_ensure_ghosting(struct mesh* m, unsigned nlayers);

#endif
