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
T* reorder_array(T const* a, unsigned const* old_to_new,
    unsigned n, unsigned width);
template <typename T>
T* reorder_array_inv(T const* a, unsigned const* new_to_old,
    unsigned n, unsigned width);
template <typename T>
T array_at(T const* a, unsigned i);
template <typename T>
void expand_into(T* out, T const* a,
    unsigned const* offsets, unsigned n, unsigned width);
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
