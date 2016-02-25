#ifndef COMPRESS_H
#define COMPRESS_H

extern unsigned const can_compress;

void* my_compress(
    void const* in_data,
    unsigned long in_size,
    unsigned long* out_size);

void* my_decompress(
    void const* in_data,
    unsigned long in_size,
    unsigned long out_size);

#endif
