#ifndef LABEL_H
#define LABEL_H

struct label {
  char* name;
  unsigned* data;
};

struct labels {
  unsigned n;
  struct label** at;
};

#endif
