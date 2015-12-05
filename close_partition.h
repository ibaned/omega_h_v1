#ifndef CLOSE_PARTITION_H
#define CLOSE_PARTITION_H

struct exchanger;

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

#endif
