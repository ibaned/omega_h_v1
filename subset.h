#ifndef SUBSET_H
#define SUBSET_H

unsigned* ints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* set_offsets);

struct subgraph {
  unsigned nverts;
  unsigned padding__;
  unsigned* offsets;
  unsigned* edges;
};

struct subgraph induce_subgraph(
    unsigned nverts,
    unsigned const* offsets,
    unsigned const* edges,
    unsigned const* set_offsets);

#endif
