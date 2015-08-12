#ifndef COLLECT_KEYS_H
#define COLLECT_KEYS_H

unsigned* collect_keys(
    unsigned elem_dim,
    unsigned key_dim,
    unsigned nkeys,
    unsigned const* elems_of_keys_offsets,
    unsigned const* elems_of_keys,
    unsigned const* elems_of_keys_directions,
    unsigned const* bad_elems,
    unsigned const* key_of_elems);

#endif
