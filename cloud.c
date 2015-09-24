#include "cloud.h"
#include "loop.h"
#include <string.h>

struct cloud {
  unsigned count;
  int padding__;
  struct tags tags;
};

struct cloud* new_cloud(unsigned count)
{
  struct cloud* c = loop_host_malloc(sizeof(*c));
  memset(c, 0, sizeof(*c));
  c->count = count;
  return c;
}

void free_cloud(struct cloud* c)
{
  free_tags(&c->tags);
}

unsigned cloud_count(struct cloud* c)
{
  return c->count;
}

struct const_tag* cloud_add_tag(struct cloud* c, enum tag_type type,
    char const* name, unsigned ncomps, void* data)
{
  return add_tag(&c->tags, type, name, ncomps, data);
}

void cloud_free_tag(struct cloud* c, char const* name)
{
  remove_tag(&c->tags, name);
}

struct const_tag* cloud_find_tag(struct cloud* c, char const* name)
{
  return (struct const_tag*) find_tag(&c->tags, name);
}

unsigned cloud_count_tags(struct cloud* c)
{
  return c->tags.n;
}

struct const_tag* cloud_get_tag(struct cloud* c, unsigned i)
{
  return get_tag(&c->tags, i);
}
