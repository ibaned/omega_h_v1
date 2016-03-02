#ifndef SHUFFLE_MESH_H
#define SHUFFLE_MESH_H

struct mesh;

unsigned* number_ents(struct mesh* m,
    unsigned ent_dim, unsigned const* vert_num);

#endif
