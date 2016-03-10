#ifndef GRAPH_HPP
#define GRAPH_HPP

struct graph {
  unsigned* offsets;
  unsigned* adj;
};

struct const_graph {
  unsigned const* const offsets;
  unsigned const* const adj;
};

struct graph* osh_new_graph(unsigned* offsets, unsigned* adj);
void osh_free_graph(struct graph* g);

#endif
