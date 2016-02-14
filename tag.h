#ifndef TAG_H
#define TAG_H

enum tag_type {
  TAG_U8,
  TAG_U32,
  TAG_U64,
  TAG_F64,
};

#define TAG_TYPES (TAG_F64+1)

struct const_tag {
  char const* const name;
  unsigned const ncomps;
  enum tag_type const type;
  union {
    void* raw;
    unsigned char* u8;
    unsigned* u32;
    unsigned long* u64;
    double* f64;
  } d;
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
void rename_tag(struct tags* ts, char const* oldname, char const* newname);
void modify_tag(struct tags* ts, char const* name, void* data);

unsigned tag_size(enum tag_type t);

void copy_tags(struct tags* a, struct tags* b, unsigned n);

struct exchanger;

void push_tag(struct exchanger* ex, struct const_tag* t, struct tags* into);
void push_tags(struct exchanger* push, struct tags* from, struct tags* into);

#endif
