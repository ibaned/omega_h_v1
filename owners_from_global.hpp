#ifndef OWNERS_FROM_GLOBAL_H
#define OWNERS_FROM_GLOBAL_H

void owners_from_global(
    unsigned n,
    unsigned long const* global_in,
    unsigned** p_own_parts,
    unsigned** p_own_idxs);

void own_idxs_from_global(
    unsigned n,
    unsigned long const* global_in,
    unsigned const* own_ranks_in,
    unsigned** p_own_idxs);

#endif
