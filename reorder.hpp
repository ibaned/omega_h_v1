#ifndef REORDER_HPP
#define REORDER_HPP

struct mesh;

unsigned* compute_ordering(struct mesh* m);

unsigned* number_ents(struct mesh* m,
    unsigned ent_dim, unsigned const* vert_num);

#endif
