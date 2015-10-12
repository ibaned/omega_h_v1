#ifndef GLOBAL_H
#define GLOBAL_H

struct mesh;

void mesh_number_simply(struct mesh* m);
unsigned long* globalize_offsets(unsigned* local, unsigned n);

#endif
