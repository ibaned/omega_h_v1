#ifndef FIELD_H
#define FIELD_H

struct field {
  char* name;
  unsigned ncomps;
  double* data;
};

struct field* new_field(char const* name, unsigned ncomps);

struct fields {
  unsigned n;
  struct field** at;
};

void add_field(struct fields* fs, struct field* f);
void free_fields(struct fields* fs);
struct field* find_field(struct fields* fs, char const* name);

#endif
