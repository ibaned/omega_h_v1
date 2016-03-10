#ifndef FILES_H
#define FILES_H

#include <stdio.h>

#include "loop.hpp"

typedef char line_t[1024];

enum endian {
  /* BIG_ENDIAN and others are often defined by
     the OS/compiler */
  MY_BIG_ENDIAN,
  MY_LITTLE_ENDIAN
};

void split_pathname(char const* pathname, char* buf,
    unsigned buf_size, char** filename, char** suffix);
char* add_enum_suffix(char const* prefix, unsigned nitems,
    unsigned item, char* buf, unsigned long buf_size);
void enum_pathname(char const* prefix, unsigned npieces,
    unsigned piece, char const* suffix, char* buf, unsigned buf_size);
void safe_scanf(FILE* f, int nitems, char const* format, ...);
void safe_read(void* p, unsigned long size, unsigned long nitems, FILE* f);
void safe_seek(FILE* f, long offset, int whence);
void seek_prefix(FILE* f,
    char* line, unsigned line_size, char const* prefix);
enum endian endianness(void);
void* generic_swap_if_needed(enum endian e, unsigned n, unsigned width,
    void const* a);
FILE* safe_fopen(char const* filename, char const* mode);

static inline LOOP_INOUT void swap_one(void* a, unsigned width)
{
  unsigned char* b = (unsigned char*) a;
  for (unsigned j = 0; j < width/2; ++j) {
    unsigned char tmp = b[j];
    b[j] = b[width - j - 1];
    b[width - j - 1] = tmp;
  }
}

#endif
