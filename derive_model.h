#ifndef DERIVE_MODEL
#define DERIVE_MODEL

struct mesh;

void mesh_derive_class_dim(struct mesh* m, double crease_angle);
void mesh_derive_class_id(struct mesh* m);
void mesh_derive_model(struct mesh* m, double crease_angle);

#endif
