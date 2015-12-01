#ifndef ARRAYS_H
#define ARRAYS_H

unsigned char* uchars_copy(unsigned char const* a, unsigned n);
unsigned* uints_copy(unsigned const* a, unsigned n);
unsigned long* ulongs_copy(unsigned long const* a, unsigned n);
double* doubles_copy(double const* a, unsigned n);

unsigned* uints_shuffle(unsigned n, unsigned const* a,
    unsigned width, unsigned const* out_of_in);
unsigned* uints_unshuffle(unsigned n, unsigned const* a,
    unsigned width, unsigned const* out_of_in);
unsigned long* ulongs_shuffle(unsigned n, unsigned long const* a,
    unsigned width, unsigned const* out_of_in);
unsigned long* ulongs_unshuffle(unsigned n, unsigned long const* a,
    unsigned width, unsigned const* out_of_in);
double* doubles_shuffle(unsigned n, double const* a,
    unsigned width, unsigned const* out_of_in);
double* doubles_unshuffle(unsigned n, double const* a,
    unsigned width, unsigned const* out_of_in);

unsigned* uints_expand(unsigned n, unsigned const* a,
    unsigned width, unsigned const* offsets);
unsigned long* ulongs_expand(unsigned n, unsigned long const* a,
    unsigned width, unsigned const* offsets);
double* doubles_expand(unsigned n, double const* a,
    unsigned width, unsigned const* offsets);

#endif
