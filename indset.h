#ifndef INDSET_H
#define INDSET_H

struct mesh;

unsigned* mesh_find_indset(struct mesh* m, unsigned ent_dim,
    unsigned const* candidates, double const* qualities);

unsigned* mesh_indset_offsets(struct mesh* m, unsigned ent_dim,
    unsigned const* candidates, double const* qualities);

#endif
