#include "field.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

struct field* new_field(char const* name, unsigned ncomps, double* data)
{
  struct field* f = malloc(sizeof(*f));
  f->name = malloc(strlen(name) + 1);
  strcpy(f->name, name);
  f->ncomps = ncomps;
  f->data = data;
  return f;
}

void add_field(struct fields* fs, struct field* f)
{
  assert(!find_field(fs, f->name));
  fs->n++;
  fs->at = realloc(fs->at, sizeof(struct field*) * fs->n);
  fs->at[fs->n - 1] = f;
}

static void free_field(struct field* f)
{
  free(f->name);
  free(f->data);
  free(f);
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
