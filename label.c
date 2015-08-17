#include "label.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

struct label* new_label(char const* name, unsigned* data)
{
  struct label* l = malloc(sizeof(*l));
  l->name = malloc(strlen(name) + 1);
  strcpy(l->name, name);
  l->data = data;
  return l;
}

void add_label(struct labels* ls, struct label* l)
{
  assert(!find_label(ls, l->name));
  ls->n++;
  ls->at = realloc(ls->at, sizeof(struct label*) * ls->n);
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
  free(l->name);
  free(l->data);
  free(l);
}

void free_labels(struct labels* ls)
{
  for (unsigned i = 0; i < ls->n; ++i)
    free_label(ls->at[i]);
  free(ls->at);
}
