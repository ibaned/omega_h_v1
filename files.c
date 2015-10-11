#include "files.h"

#include <assert.h>
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

void split_filename(char const* filename, char* buf,
    unsigned buf_size, char** suffix)
{
  assert(strlen(filename) < buf_size);
  strcpy(buf, filename);
  char* dot = strrchr(buf, '.');
  assert(dot);
  *dot = '\0';
  if (suffix)
    *suffix = dot + 1;
}

void parallel_filename(char const* prefix, unsigned npieces,
    unsigned piece, char const* suffix, char* buf, unsigned buf_size)
{
  unsigned ndig = count_digits(npieces);
  unsigned long prelen = strlen(prefix);
  unsigned long suflen = strlen(suffix);
  assert(prelen + ndig + suflen < buf_size);
  memcpy(buf, prefix, prelen);
  buf += prelen;
  if (ndig) {
    sprintf(buf, "%0*u", (int) ndig, piece);
    buf += ndig;
  }
  memcpy(buf, suffix, suflen);
}
