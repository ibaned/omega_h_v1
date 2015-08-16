#ifndef GRAPH_H
#define GRAPH_H

struct graph {
  unsigned* offsets;
  unsigned* adj;
};

struct const_graph {
  unsigned const* const offsets;
  unsigned const* const adj;
};

struct graph* new_graph(unsigned* offsets, unsigned* adj);
void free_graph(struct graph* g);

#endif
