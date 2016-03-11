#include "tag.hpp"

#include <cassert>
#include <cstring>

#include "arrays.hpp"
#include "exchanger.hpp"
#include "loop.hpp"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

struct tag {
  char* name;
  unsigned ncomps;
  enum tag_type type;
  enum osh_transfer transfer_type;
  void* data;
};

#ifdef __clang__
#pragma clang diagnostic pop
#endif

unsigned tag_size(enum tag_type t)
{
  switch (t) {
    case TAG_U8:  return sizeof(unsigned char);
    case TAG_U32: return sizeof(unsigned);
    case TAG_U64: return sizeof(unsigned long);
    case TAG_F64: return sizeof(double);
  }
  LOOP_NORETURN(0);
}

static struct tag* new_tag(char const* name, enum tag_type type,
    unsigned ncomps, enum osh_transfer tt, void* data)
{
  struct tag* t = LOOP_HOST_MALLOC(struct tag, 1);
  t->name = LOOP_HOST_MALLOC(char, strlen(name) + 1);
  strcpy(t->name, name);
  t->ncomps = ncomps;
  t->type = type;
  t->transfer_type = tt;
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
  for (unsigned i = 0; i < ts->n; ++i)
    free_tag(ts->at[i]);
  loop_host_free(ts->at);
}

struct const_tag* add_tag(struct tags* ts, enum tag_type type, char const* name,
    unsigned ncomps, void* data)
{
  return add_tag2(ts, type, name, ncomps, OSH_TRANSFER_NOT, data);
}

struct const_tag* add_tag2(struct tags* ts, enum tag_type type, char const* name,
    unsigned ncomps, enum osh_transfer tt, void* data)
{
  assert(!find_tag(ts, name));
  struct tag* t = new_tag(name, type, ncomps, tt, data);
  if (ts->n == ts->cap) {
    ts->cap = (3 * ts->cap) / 2;
    if (ts->cap < 10)
      ts->cap = 10;
    ts->at = LOOP_HOST_REALLOC(struct tag*, ts->at, ts->cap);
  }
  ts->at[ts->n++] = t;
  return reinterpret_cast<struct const_tag*>(t);
}

static unsigned find_i(struct tags* ts, char const* name)
{
  unsigned i;
  for (i = 0; i < ts->n; ++i)
    if (!strcmp(ts->at[i]->name, name))
      break;
  return i;
}

void remove_tag(struct tags* ts, char const* name)
{
  for (unsigned type = 0; type < TAG_TYPES; ++type) {
    unsigned i = find_i(ts, name);
    if (i >= ts->n)
      continue;
    free_tag(ts->at[i]);
    ts->n--;
    for (; i < ts->n; ++i)
      ts->at[i] = ts->at[i + 1];
    return;
  }
  assert(0);
}

struct const_tag* find_tag(struct tags* ts, char const* name)
{
  unsigned i = find_i(ts, name);
  if (i < ts->n)
    return get_tag(ts, i);
  return 0;
}

unsigned count_tags(struct tags* ts)
{
  return ts->n;
}

struct const_tag* get_tag(struct tags* ts, unsigned i)
{
  return reinterpret_cast<struct const_tag*>(ts->at[i]);
}

void rename_tag(struct tags* ts, char const* oldname, char const* newname)
{
  if (!strcmp(oldname, newname))
    return;
  assert (!find_tag(ts, newname));
  unsigned i = find_i(ts, oldname);
  assert(i < ts->n);
  loop_host_free(ts->at[i]->name);
  ts->at[i]->name = LOOP_HOST_MALLOC(char, strlen(newname) + 1);
  strcpy(ts->at[i]->name, newname);
}

void modify_tag(struct tags* ts, char const* name, void* data)
{
  unsigned i = find_i(ts, name);
  struct tag* t = ts->at[i];
  loop_free(t->data);
  t->data = data;
}

void copy_tags(struct tags* a, struct tags* b, unsigned n)
{
  for (unsigned i = 0; i < count_tags(a); ++i) {
    struct const_tag* t = get_tag(a, i);
    void* data = 0;
    switch (t->type) {
      case TAG_U8:  data = uchars_copy(t->d.u8, n * t->ncomps);
                    break;
      case TAG_U32: data = uints_copy(t->d.u32, n * t->ncomps);
                    break;
      case TAG_U64: data = ulongs_copy(t->d.u64, n * t->ncomps);
                    break;
      case TAG_F64: data = doubles_copy(t->d.f64, n * t->ncomps);
                    break;
    };
    add_tag2(b, t->type, t->name, t->ncomps, t->transfer_type, data);
  }
}

void push_tag(struct exchanger* ex, struct const_tag* t, struct tags* into)
{
  void* data_out = 0;
  switch (t->type) {
    case TAG_U8:
      break;
    case TAG_U32:
      data_out = exchange_uints(ex, t->ncomps, t->d.u32, EX_FOR, EX_ROOT);
      break;
    case TAG_U64:
      data_out = exchange_ulongs(ex, t->ncomps, t->d.u64, EX_FOR, EX_ROOT);
      break;
    case TAG_F64:
      data_out = exchange_doubles(ex, t->ncomps, t->d.f64, EX_FOR, EX_ROOT);
      break;
  }
  if (find_tag(into, t->name))
    modify_tag(into, t->name, data_out);
  else
    add_tag2(into, t->type, t->name, t->ncomps, t->transfer_type, data_out);
}

void push_tags(struct exchanger* ex, struct tags* from, struct tags* into)
{
  for (unsigned i = 0; i < count_tags(from); ++i)
    push_tag(ex, get_tag(from, i), into);
}
