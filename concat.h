#ifndef CONCAT_H
#define CONCAT_H

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
