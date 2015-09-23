#ifndef FIELD_H
#define FIELD_H

struct const_field {
  char const* const name;
  unsigned const ncomps;
  int padding__;
  double const* const data;
};

struct fields {
  unsigned n;
  int padding__;
  struct field** at;
};

struct const_field* add_field(struct fields* fs, char const* name,
    unsigned ncomps, double* data);
void remove_field(struct fields* fs, char const* name);
void free_fields(struct fields* fs);
struct const_field* find_field(struct fields* fs, char const* name);
struct const_field* get_field(struct fields* fs, unsigned i);

#endif
