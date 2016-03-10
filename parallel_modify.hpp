#ifndef PARALLEL_MODIFY
#define PARALLEL_MODIFY

struct mesh;

void set_own_ranks_by_indset(
    struct mesh* m,
    unsigned key_dim);

void inherit_globals(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned ent_dim,
    unsigned const* offset_of_same_ents);

#endif
