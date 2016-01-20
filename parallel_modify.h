#ifndef PARALLEL_MODIFY
#define PARALLEL_MODIFY

struct mesh;

void set_own_ranks_by_indset(
    struct mesh* m,
    unsigned key_dim);

#endif
