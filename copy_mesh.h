#ifndef COPY_MESH_H
#define COPY_MESH_H

struct tags;
struct mesh;

void copy_tags(struct tags* a, struct tags* b, unsigned n);
struct mesh* copy_mesh(struct mesh* a);

#endif
