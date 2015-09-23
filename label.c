#include "label.h"
#include "loop.h"
#include <string.h>
#include <assert.h>

struct label {
  char* name;
  unsigned* data;
};

static struct label* new_label(char const* name, unsigned* data)
{
  struct label* l = loop_host_malloc(sizeof(*l));
  l->name = loop_host_malloc(strlen(name) + 1);
  strcpy(l->name, name);
  l->data = data;
  return l;
}

static void free_label(struct label* l)
{
  loop_host_free(l->name);
  loop_free(l->data);
  loop_host_free(l);
}

struct const_label* add_label(struct labels* ls, char const* name,
    unsigned* data)
{
  assert(!find_label(ls, name));
  struct label* l = new_label(name, data);
  ls->n++;
  ls->at = loop_host_realloc(ls->at, sizeof(struct label*) * ls->n);
  ls->at[ls->n - 1] = l;
  return (struct const_label*) l;
}

static unsigned find_label_i(struct labels* ls, char const* name)
{
  unsigned i;
  for (i = 0; i < ls->n; ++i)
    if (!strcmp(ls->at[i]->name, name))
      break;
  return i;
}

void remove_label(struct labels* ls, char const* name)
{
  unsigned i = find_label_i(ls, name);
  assert(i < ls->n);
  free_label(ls->at[i]);
  ls->n--;
  for (; i < ls->n; ++i)
    ls->at[i] = ls->at[i + 1];
}

void free_labels(struct labels* ls)
{
  for (unsigned i = 0; i < ls->n; ++i)
    free_label(ls->at[i]);
  loop_host_free(ls->at);
}

struct const_label* find_label(struct labels* ls, char const* name)
{
  unsigned i = find_label_i(ls, name);
  if (i < ls->n)
    return get_label(ls, (unsigned) i);
  return 0;
}

struct const_label* get_label(struct labels* ls, unsigned i)
{
  return (struct const_label*) ls->at[i];
}
