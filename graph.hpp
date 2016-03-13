#ifndef GRAPH_HPP
#define GRAPH_HPP

namespace omega_h {

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

}

#endif
