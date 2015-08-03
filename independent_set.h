#ifndef INDEPENDENT_SET_H
#define INDEPENDENT_SET_H

unsigned* find_independent_set(
    unsigned nverts,
    unsigned const* offsets,
    unsigned const* edges,
    unsigned const* filter,
    double const* goodness);

#endif
