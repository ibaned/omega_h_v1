#include "field.h"
#include <assert.h>  // for assert
#include <string.h>  // for strcpy
#include "loop.h"    // for free, malloc, realloc
#include <string.h>  // for strcmp, strlen

struct field* new_field(char const* name, unsigned ncomps, double* data)
{
  struct field* f = loop_malloc(sizeof(*f));
  f->name = loop_host_malloc(strlen(name) + 1);
  strcpy(f->name, name);
  f->ncomps = ncomps;
  f->data = data;
  return f;
}

void free_field(struct field* f)
{
  loop_host_free(f->name);
  loop_free(f->data);
  loop_host_free(f);
}

void add_field(struct fields* fs, struct field* f)
{
  assert(!find_field(fs, f->name));
  fs->n++;
  fs->at = loop_host_realloc(fs->at, sizeof(struct field*) * fs->n);
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
  loop_host_free(fs->at);
}

struct field* find_field(struct fields* fs, char const* name)
{
  for (unsigned i = 0; i < fs->n; ++i)
    if (!strcmp(fs->at[i]->name, name))
      return fs->at[i];
  return 0;
}
