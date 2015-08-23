#ifndef FIELD_H
#define FIELD_H

struct field {
  char* name;
  unsigned ncomps;
  double* data;
};

struct const_field {
  char const* const name;
  unsigned const ncomps;
  double const* const data;
};

struct field* new_field(char const* name, unsigned ncomps, double* data);
void free_field(struct field* f);

struct fields {
  unsigned n;
  struct field** at;
};

void add_field(struct fields* fs, struct field* f);
void remove_field(struct fields* fs, struct field* f);
void free_fields(struct fields* fs);
struct field* find_field(struct fields* fs, char const* name);

#endif
