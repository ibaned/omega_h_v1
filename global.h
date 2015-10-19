#ifndef GLOBAL_H
#define GLOBAL_H

#define MAX_NEIGHBORS 128

struct mesh;

void mesh_number_simply(struct mesh* m);
unsigned long* globalize_offsets(unsigned* local, unsigned n);
void global_to_linpart(unsigned long* global, unsigned n,
    unsigned long total, unsigned nparts,
    unsigned** p_part, unsigned** p_local);

#endif
