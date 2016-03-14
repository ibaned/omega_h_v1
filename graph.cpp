#include "graph.hpp"

#include "loop.hpp"

namespace omega_h {

struct graph* new_graph(unsigned* offsets, unsigned* adj)
{
  struct graph* g = LOOP_HOST_MALLOC(struct graph, 1);
  g->offsets = offsets;
  g->adj = adj;
  return g;
}

void free_graph(struct graph* g)
{
  if (!g)
    return;
  loop_free(g->offsets);
  loop_free(g->adj);
  loop_host_free(g);
}

}
