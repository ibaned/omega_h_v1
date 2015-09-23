#include "field.h"
#include <assert.h>  // for assert
#include <string.h>  // for strcpy
#include "loop.h"    // for free, malloc, realloc
#include <string.h>  // for strcmp, strlen

struct field {
  char* name;
  unsigned ncomps;
  int padding__;
  double* data;
};

static struct field* new_field(char const* name, unsigned ncomps, double* data)
{
  struct field* f = loop_host_malloc(sizeof(*f));
  f->name = loop_host_malloc(strlen(name) + 1);
  strcpy(f->name, name);
  f->ncomps = ncomps;
  f->data = data;
  return f;
}

static void free_field(struct field* f)
{
  loop_host_free(f->name);
  loop_free(f->data);
  loop_host_free(f);
}

struct const_field* add_field(struct fields* fs, char const* name,
    unsigned ncomps, double* data)
{
  assert(!find_field(fs, name));
  struct field* f = new_field(name, ncomps, data);
  fs->n++;
  fs->at = loop_host_realloc(fs->at, sizeof(struct field*) * fs->n);
  fs->at[fs->n - 1] = f;
  return (struct const_field*) f;
}

static unsigned find_field_i(struct fields* fs, char const* name)
{
  unsigned i;
  for (i = 0; i < fs->n; ++i)
    if (!strcmp(fs->at[i]->name, name))
      break;
  return i;
}

void remove_field(struct fields* fs, char const* name)
{
  unsigned i = find_field_i(fs, name);
  assert(i < fs->n);
  free_field(fs->at[i]);
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

struct const_field* find_field(struct fields* fs, char const* name)
{
  unsigned i = find_field_i(fs, name);
  if (i < fs->n)
    return get_field(fs, (unsigned) i);
  return 0;
}

struct const_field* get_field(struct fields* fs, unsigned i)
{
  return (struct const_field*) fs->at[i];
}
