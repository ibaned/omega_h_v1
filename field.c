#include "field.h"
#include <assert.h>  // for assert
#include <string.h>  // for strcpy
#include <stdlib.h>  // for free, malloc, realloc
#include <string.h>  // for strcmp, strlen

struct field* new_field(char const* name, unsigned ncomps, double* data)
{
  struct field* f = malloc(sizeof(*f));
  f->name = malloc(strlen(name) + 1);
  strcpy(f->name, name);
  f->ncomps = ncomps;
  f->data = data;
  return f;
}

void free_field(struct field* f)
{
  free(f->name);
  free(f->data);
  free(f);
}

void add_field(struct fields* fs, struct field* f)
{
  assert(!find_field(fs, f->name));
  fs->n++;
  fs->at = realloc(fs->at, sizeof(struct field*) * fs->n);
  fs->at[fs->n - 1] = f;
}

void remove_field(struct fields* fs, struct field* f)
{
  unsigned i;
  for (i = 0; i < fs->n; ++i)
    if (fs->at[i] == f)
      break;
  fs->n--;
  for (; i < fs->n; ++i)
    fs->at[i] = fs->at[i + 1];
}

void free_fields(struct fields* fs)
{
  for (unsigned i = 0; i < fs->n; ++i)
    free_field(fs->at[i]);
  free(fs->at);
}

struct field* find_field(struct fields* fs, char const* name)
{
  for (unsigned i = 0; i < fs->n; ++i)
    if (!strcmp(fs->at[i]->name, name))
      return fs->at[i];
  return 0;
}
