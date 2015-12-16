#ifndef CLOSE_PARTITION_H
#define CLOSE_PARTITION_H

struct exchanger;
struct mesh;

void close_partition_exchangers(
    unsigned nacopies,
    unsigned nbowners,
    unsigned const* buses_by_acopies_offsets,
    unsigned const* buse_own_ranks,
    unsigned const* buse_own_ids,
    struct exchanger** p_buse_to_own,
    struct exchanger** p_bown_to_copy);

void close_partition(
    unsigned nacopies,
    unsigned nbowners,
    unsigned const* buses_by_acopies_offsets,
    unsigned const* buse_own_ranks,
    unsigned const* buse_own_ids,
    unsigned* p_nbcopies,
    unsigned** p_bcopy_own_ranks,
    unsigned** p_bcopy_own_ids);

void get_down_use_owners(
    struct mesh* m,
    unsigned high_dim,
    unsigned low_dim,
    unsigned** p_use_own_ranks,
    unsigned** p_use_own_ids,
    unsigned** p_uses_of_highs_offsets);

void push_use_owners(
    struct exchanger* push,
    unsigned const* use_own_ranks_in,
    unsigned const* use_own_ids_in,
    unsigned const* offsets_in,
    unsigned** p_use_own_ranks_out,
    unsigned** p_use_own_ids_out,
    unsigned** p_offsets_out);

#endif
