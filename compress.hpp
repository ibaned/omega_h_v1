#ifndef COMPRESS_HPP
#define COMPRESS_HPP

namespace omega_h {

extern unsigned const can_compress;

void* my_compress(
    void const* in_data,
    unsigned long in_size,
    unsigned long* out_size);

void* my_decompress(
    void const* in_data,
    unsigned long in_size,
    unsigned long out_size);

}

#endif
