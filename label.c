#include "label.h"
#include "loop.h"
#include <string.h>
#include <assert.h>

struct label* new_label(char const* name, unsigned* data)
{
  struct label* l = loop_malloc(sizeof(*l));
  l->name = loop_host_malloc(strlen(name) + 1);
  strcpy(l->name, name);
  l->data = data;
  return l;
}

void add_label(struct labels* ls, struct label* l)
{
  assert(!find_label(ls, l->name));
  ls->n++;
  ls->at = loop_host_realloc(ls->at, sizeof(struct label*) * ls->n);
  ls->at[ls->n - 1] = l;
}

struct label* find_label(struct labels* ls, char const* name)
{
  for (unsigned i = 0; i < ls->n; ++i)
    if (!strcmp(ls->at[i]->name, name))
      return ls->at[i];
  return 0;
}

static void free_label(struct label* l)
{
  loop_host_free(l->name);
  loop_free(l->data);
  loop_host_free(l);
}

void free_labels(struct labels* ls)
{
  for (unsigned i = 0; i < ls->n; ++i)
    free_label(ls->at[i]);
  loop_host_free(ls->at);
}
