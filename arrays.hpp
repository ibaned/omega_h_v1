#ifndef ARRAYS_H
#define ARRAYS_H

void uchars_memcpy(unsigned char* dst,
    unsigned char const* src, unsigned n);
void uints_memcpy(unsigned* dst,
    unsigned const* src, unsigned n);
void ulongs_memcpy(unsigned long* dst,
    unsigned long const* src, unsigned n);
void doubles_memcpy(double* dst,
    double const* src, unsigned n);
void uintptrs_memcpy(unsigned** dst,
    unsigned* const* src, unsigned n);

unsigned char* uchars_copy(unsigned char const* a, unsigned n);
unsigned* uints_copy(unsigned const* a, unsigned n);
unsigned long* ulongs_copy(unsigned long const* a, unsigned n);
double* doubles_copy(double const* a, unsigned n);
unsigned** uintptrs_copy(unsigned* const* a, unsigned n);

unsigned char* uchars_to_device(unsigned char const* a, unsigned n);
unsigned* uints_to_device(unsigned const* a, unsigned n);
unsigned long* ulongs_to_device(unsigned long const* a, unsigned n);
double* doubles_to_device(double const* a, unsigned n);
unsigned** uintptrs_to_device(unsigned* const* a, unsigned n);

unsigned char* uchars_to_host(unsigned char const* a, unsigned n);
unsigned* uints_to_host(unsigned const* a, unsigned n);
unsigned long* ulongs_to_host(unsigned long const* a, unsigned n);
double* doubles_to_host(double const* a, unsigned n);
unsigned** uintptrs_to_host(unsigned* const* a, unsigned n);

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

void uchars_expand_into(unsigned n, unsigned width,
    unsigned char const* a, unsigned const* offsets, unsigned char* out);
void uints_expand_into(unsigned n, unsigned width,
    unsigned const* a, unsigned const* offsets, unsigned* out);
void ulongs_expand_into(unsigned n, unsigned width,
    unsigned long const* a, unsigned const* offsets, unsigned long* out);
void doubles_expand_into(unsigned n, unsigned width,
    double const* a, unsigned const* offsets, double* out);

unsigned char* uchars_expand(unsigned n, unsigned width,
    unsigned char const* a, unsigned const* offsets);
unsigned* uints_expand(unsigned n, unsigned width,
    unsigned const* a, unsigned const* offsets);
unsigned long* ulongs_expand(unsigned n, unsigned width,
    unsigned long const* a, unsigned const* offsets);
double* doubles_expand(unsigned n, unsigned width,
    double const* a, unsigned const* offsets);

unsigned* concat_uints(unsigned width,
    unsigned const* a, unsigned na,
    unsigned const* b, unsigned nb);
double* concat_doubles(unsigned width,
    double const* a, unsigned na,
    double const* b, unsigned nb);

unsigned char* uchars_filled(unsigned n, unsigned char v);
unsigned* uints_filled(unsigned n, unsigned v);
unsigned long* ulongs_filled(unsigned n, unsigned long v);
double* doubles_filled(unsigned n, double v);

unsigned char uchars_at(unsigned char const* a, unsigned i);
unsigned uints_at(unsigned const* a, unsigned i);

void doubles_max_into(unsigned n, unsigned width,
    double const* a, unsigned const* offsets, double* out);

#endif
