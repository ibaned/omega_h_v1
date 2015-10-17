#include "files.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

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
