#ifndef CONCAT_H
#define CONCAT_H

unsigned* ints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* offsets);

double* doubles_subset(
    unsigned n,
    unsigned width,
    double const* a,
    unsigned const* offsets);

unsigned* concat_ints(
    unsigned narrs,
    unsigned width,
    unsigned const* sizes,
    unsigned const* const* arrs);

double* concat_doubles(
    unsigned narrs,
    unsigned width,
    unsigned const* sizes,
    double const* const* arrs);

#endif
