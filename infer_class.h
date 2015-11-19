#ifndef INFER_CLASS_H
#define INFER_CLASS_H

void infer_class(
    unsigned ent_dim,
    unsigned nents,
    unsigned const* verts_of_ents,
    unsigned const* class_dim_of_verts,
    unsigned const* class_id_of_verts,
    unsigned** p_class_dim_of_ents,
    unsigned** p_class_id_of_ents);

struct mesh;

unsigned const* mesh_ask_class_dim(struct mesh* m, unsigned dim);
unsigned const* mesh_ask_class_id(struct mesh* m, unsigned dim);

#endif
