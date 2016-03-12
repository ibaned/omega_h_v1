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
T* reorder_array(unsigned n, T const* a,
    unsigned width, unsigned const* old_to_new);
template <typename T>
T* reorder_array_inv(unsigned n, T const* a,
    unsigned width, unsigned const* new_to_old);
template <typename T>
T array_at(T const* a, unsigned i);
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
template <typename T>
T* filled_array(unsigned n, T v);
void doubles_max_into(unsigned n, unsigned width,
    double const* a, unsigned const* offsets, double* out);

#endif
