#include "compress.hpp"

#include <cassert>

#include "loop_host.hpp"

#if USE_ZLIB

#include <zlib.h>

unsigned const can_compress = 1;

void* my_compress(
    void const* in_data,
    unsigned long in_size,
    unsigned long* out_size)
{
  uLongf bufsize = compressBound(in_size);
  void* out_data = LOOP_HOST_MALLOC(unsigned char, bufsize);
  int ret = compress(
      static_cast<Bytef*>(out_data),
      &bufsize,
      static_cast<Bytef const*>(in_data),
      in_size);
  assert(ret == Z_OK);
  *out_size = bufsize;
  return out_data;
}

void* my_decompress(
    void const* in_data,
    unsigned long in_size,
    unsigned long out_size)
{
  uLongf bufsize = out_size;
  void* out_data = LOOP_HOST_MALLOC(unsigned char, out_size);
  int ret = uncompress(
      static_cast<Bytef*>(out_data),
      &bufsize,
      static_cast<Bytef const*>(in_data),
      in_size);
  assert(ret == Z_OK);
  assert(bufsize == out_size);
  return out_data;
}

#else

unsigned const can_compress = 0;

void* my_compress(
    void const* in_data,
    unsigned long in_size,
    unsigned long* out_size)
{
  *out_size = in_size;
  return loop_host_copy(in_data, in_size);
}

void* my_decompress(
    void const* in_data,
    unsigned long in_size,
    unsigned long out_size)
{
  assert(in_size == out_size);
  return loop_host_copy(in_data, in_size);
}

#endif
