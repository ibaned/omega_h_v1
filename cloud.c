#include "cloud.h"
#include "field.h"
#include "label.h"
#include "loop.h"
#include <string.h>

struct cloud {
  unsigned count;
  int padding__;
  struct fields fields;
  struct labels labels;
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
  free_fields(&c->fields);
  free_labels(&c->labels);
}

unsigned cloud_count(struct cloud* c)
{
  return c->count;
}

struct const_field* cloud_add_field(struct cloud* c, char const* name,
    unsigned ncomps, double* data)
{
  return add_field(&c->fields, name, ncomps, data);
}

void cloud_free_field(struct cloud* c, char const* name)
{
  remove_field(&c->fields, name);
}

struct const_field* cloud_find_field(struct cloud* c, char const* name)
{
  return (struct const_field*) find_field(&c->fields, name);
}

unsigned cloud_count_fields(struct cloud* c)
{
  return c->fields.n;
}

struct const_field* cloud_get_field(struct cloud* c, unsigned i)
{
  return get_field(&c->fields, i);
}

struct const_label* cloud_add_label(struct cloud* c, char const* name,
    unsigned* data)
{
  return add_label(&c->labels, name, data);
}

void cloud_free_label(struct cloud* c, char const* name)
{
  remove_label(&c->labels, name);
}

struct const_label* cloud_find_label(struct cloud* c, char const* name)
{
  return find_label(&c->labels, name);
}

unsigned cloud_count_labels(struct cloud* c)
{
  return c->labels.n;
}

struct const_label* cloud_get_label(struct cloud* c, unsigned i)
{
  return get_label(&c->labels, i);
}
