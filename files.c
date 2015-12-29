#include "files.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "loop.h"

static unsigned count_digits(unsigned x)
{
  unsigned l = 0;
  while (x) {
    ++l;
    x /= 10;
  }
  return l;
}

void split_pathname(char const* pathname, char* buf,
    unsigned buf_size, char** filename, char** suffix)
{
  assert(strlen(pathname) < buf_size);
  strcpy(buf, pathname);
  char* dot = strrchr(buf, '.');
  assert(dot);
  *dot = '\0';
  if (suffix)
    *suffix = dot + 1;
  if (!filename)
    return;
  char* slash = strrchr(buf, '/');
  if (slash) {
    *slash = '\0';
    *filename = slash + 1;
  } else {
    *filename = buf;
  }
}

char* add_enum_suffix(char const* prefix, unsigned nitems,
    unsigned item, char* buf, unsigned long buf_size)
{
  unsigned ndig = count_digits(nitems);
  unsigned long prelen = strlen(prefix);
  assert(prelen + 1 + ndig < buf_size);
  memcpy(buf, prefix, prelen);
  buf += prelen;
  if (nitems > 1) {
    sprintf(buf, "_%0*u", (int) ndig, item);
    buf += 1 + ndig;
  }
  *buf = '\0';
  return buf;
}

void enum_pathname(char const* prefix, unsigned npieces,
    unsigned piece, char const* suffix, char* buf, unsigned buf_size)
{
  unsigned long suflen = strlen(suffix);
  assert(suflen + 1 < buf_size);
  buf = add_enum_suffix(prefix, npieces, piece, buf, buf_size - suflen - 1);
  *buf = '.';
  ++buf;
  memcpy(buf, suffix, suflen);
  buf += suflen;
  *buf = '\0';
}

void safe_scanf(FILE* f, int nitems, char const* format, ...)
{
  va_list ap;
  va_start(ap, format);
  int r = vfscanf(f, format, ap);
  va_end(ap);
  assert(r == nitems);
}

void safe_read(void* p, unsigned long size, unsigned long nitems, FILE* f)
{
  unsigned long r = fread(p, size, nitems, f);
  assert(r == nitems);
}

void safe_seek(FILE* f, long offset, int whence)
{
  int ret = fseek(f, offset, whence);
  assert(ret == 0);
}

void seek_prefix(FILE* f,
    char* line, unsigned line_size, char const* prefix)
{
  unsigned long pl = strlen(prefix);
  while (fgets(line, (int) line_size, f))
    if (!strncmp(line, prefix, pl))
      return;
  assert(0);
}

enum endian endianness(void)
{
  static unsigned short canary = 0x1;
  unsigned char* p = (unsigned char*) (&canary);
  if (*p == 0x1)
    return MY_LITTLE_ENDIAN;
  return MY_BIG_ENDIAN;
}

LOOP_KERNEL(swap_kern, unsigned width, unsigned char* b)
  for (unsigned j = 0; j < width/2; ++j) {
    unsigned char tmp = b[i * width + j];
    b[i * width + j] = b[i * width + (width - j - 1)];
    b[i * width + (width - j - 1)] = tmp;
  }
}

void* generic_swap_if_needed(enum endian e, unsigned n, unsigned width,
    void const* a)
{
  unsigned char* b = loop_to_device(a, n * width);
  if (e != endianness() && width > 1) {
    assert(width % 2 == 0);
    LOOP_EXEC(swap_kern, n, width, b);
  }
  return b;
}
