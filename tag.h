#ifndef TAG_H
#define TAG_H

enum tag_type {
  U32,
  F64,
  TAG_TYPES
};

struct const_tag {
  char const* const name;
  unsigned const ncomps;
  int const padding__;
  void const* const data;
};

struct type_tags {
  unsigned n;
  int padding__;
  struct tag** at;
};

struct tags {
  struct type_tags of[TAG_TYPES];
};

void free_tags(struct tags* ts);
struct const_tag* add_tag(struct tags* ts, enum tag_type type, char const* name,
    unsigned ncomps, void* data);
void remove_tag(struct tags* ts, enum tag_type type, char const* name);
struct const_tag* find_tag(struct tags* ts, char const* name);
unsigned count_tags(struct tags* ts, enum tag_type type);
struct const_tag* get_tag(struct tags* ts, enum tag_type type, unsigned i);

#endif
