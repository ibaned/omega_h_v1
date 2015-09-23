#include "tag.h"
#include <assert.h>
#include <string.h>
#include "loop.h"

struct tag {
  char* name;
  unsigned ncomps;
  int padding__;
  void* data;
};

static struct tag* new_tag(char const* name, unsigned ncomps, void* data)
{
  struct tag* t = loop_host_malloc(sizeof(*t));
  t->name = loop_host_malloc(strlen(name) + 1);
  strcpy(t->name, name);
  t->ncomps = ncomps;
  t->data = data;
  return t;
}

static void free_tag(struct tag* t)
{
  loop_host_free(t->name);
  loop_free(t->data);
  loop_host_free(t);
}

void free_tags(struct tags* ts)
{
  for (enum tag_type type = 0; type < TAG_TYPES; ++type) {
    for (unsigned i = 0; i < ts->of[type].n; ++i)
      free_tag(ts->of[type].at[i]);
    loop_host_free(ts->of[type].at);
  }
}

struct const_tag* add_tag(struct tags* ts, enum tag_type type, char const* name,
    unsigned ncomps, void* data)
{
  assert(!find_tag(ts, name));
  struct tag* t = new_tag(name, ncomps, data);
  struct type_tags* tt = &ts->of[type];
  tt->n++;
  tt->at = loop_host_realloc(tt->at, sizeof(struct tag*) * tt->n);
  tt->at[tt->n - 1] = t;
  return (struct const_tag*) t;
}

static unsigned find_i(struct type_tags* tt, char const* name)
{
  unsigned i;
  for (i = 0; i < tt->n; ++i)
    if (!strcmp(tt->at[i]->name, name))
      break;
  return i;
}

void remove_tag(struct tags* ts, enum tag_type type, char const* name)
{
  struct type_tags* tt = &ts->of[type];
  unsigned i = find_i(tt, name);
  assert(i < tt->n);
  free_tag(tt->at[i]);
  tt->n--;
  for (; i < tt->n; ++i)
    tt->at[i] = tt->at[i + 1];
}

struct const_tag* find_tag(struct tags* ts, char const* name)
{
  for (enum tag_type type = 0; type < TAG_TYPES; ++type) {
    struct type_tags* tt = &ts->of[type];
    unsigned i = find_i(tt, name);
    if (i < tt->n)
      return get_tag(ts, type, i);
  }
  return 0;
}

unsigned count_tags(struct tags* ts, enum tag_type type)
{
  return ts->of[type].n;
}

struct const_tag* get_tag(struct tags* ts, enum tag_type type, unsigned i)
{
  return (struct const_tag*) ts->of[type].at[i];
}
