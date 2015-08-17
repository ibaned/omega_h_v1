#ifndef LABEL_H
#define LABEL_H

struct label {
  char* name;
  unsigned* data;
};

struct const_label {
  char const* const name;
  unsigned const* const data;
};

struct labels {
  unsigned n;
  struct label** at;
};

struct label* new_label(char const* name, unsigned* data);
void add_label(struct labels* ls, struct label* l);
struct label* find_label(struct labels* ls, char const* name);
void free_labels(struct labels* ls);

#endif
