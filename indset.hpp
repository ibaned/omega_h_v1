#ifndef INDSET_HPP
#define INDSET_HPP

struct mesh;

unsigned* mesh_find_indset(struct mesh* m, unsigned ent_dim,
    unsigned const* candidates, double const* qualities);

unsigned* mesh_indset_offsets(struct mesh* m, unsigned ent_dim,
    unsigned const* candidates, double const* qualities);

#endif
