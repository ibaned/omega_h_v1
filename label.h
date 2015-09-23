#ifndef LABEL_H
#define LABEL_H

struct const_label {
  char const* const name;
  unsigned const* const data;
};

struct labels {
  unsigned n;
  int padding_;
  struct label** at;
};

struct const_label* add_label(struct labels* ls, char const* name,
    unsigned* data);
void remove_label(struct labels* ls, char const* name);
void free_labels(struct labels* ls);
struct const_label* find_label(struct labels* ls, char const* name);
struct const_label* get_label(struct labels* ls, unsigned i);

#endif
