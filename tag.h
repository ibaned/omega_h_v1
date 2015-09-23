#ifndef TAG_H
#define TAG_H

enum tag_type {
  TAG_U32,
  TAG_F64,
};

#define TAG_TYPES (TAG_F64+1)

struct const_tag {
  char const* const name;
  unsigned const ncomps;
  enum tag_type const type;
  void const* const data;
};

struct tags {
  unsigned n;
  int padding__;
  struct tag** at;
};

void free_tags(struct tags* ts);
struct const_tag* add_tag(struct tags* ts, enum tag_type type, char const* name,
    unsigned ncomps, void* data);
void remove_tag(struct tags* ts, char const* name);
struct const_tag* find_tag(struct tags* ts, char const* name);
unsigned count_tags(struct tags* ts);
struct const_tag* get_tag(struct tags* ts, unsigned i);

#endif
