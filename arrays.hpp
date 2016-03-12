#ifndef ARRAYS_HPP
#define ARRAYS_HPP

template <typename T>
void array_memcpy(T* dst, T const* src, unsigned n);

template <typename T>
T* copy_array(T const* a, unsigned n);

template <typename T>
T* array_to_device(T const* a, unsigned n);

template <typename T>
T* array_to_host(T const* a, unsigned n);

template <typename T>
T* shuffle_array(unsigned n, T const* a,
    unsigned width, unsigned const* out_of_in);
template <typename T>
T* unshuffle_array(unsigned n, T const* a,
    unsigned width, unsigned const* out_of_in);

template <typename T>
void expand_into(unsigned n, unsigned width,
    T const* a, unsigned const* offsets, T* out);
template <typename T>
T* expand_array(unsigned n, unsigned width,
    T const* a, unsigned const* offsets);

template <typename T>
T* concat_arrays(unsigned width,
    T const* a, unsigned na,
    T const* b, unsigned nb);

unsigned char* uchars_filled(unsigned n, unsigned char v);
unsigned* uints_filled(unsigned n, unsigned v);
unsigned long* ulongs_filled(unsigned n, unsigned long v);
double* doubles_filled(unsigned n, double v);

unsigned char uchars_at(unsigned char const* a, unsigned i);
unsigned uints_at(unsigned const* a, unsigned i);

void doubles_max_into(unsigned n, unsigned width,
    double const* a, unsigned const* offsets, double* out);

#endif
