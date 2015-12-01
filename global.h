#ifndef GLOBAL_H
#define GLOBAL_H

struct mesh;

void mesh_number_simply(struct mesh* m);
unsigned long* globalize_offsets(unsigned* local, unsigned n);
void global_to_linpart(unsigned long const* global, unsigned n,
    unsigned long total, unsigned nparts,
    unsigned** p_part, unsigned** p_local);
unsigned linpart_size(unsigned long total, unsigned nparts, unsigned part);

#endif
