#include "graph.h"

#include "loop.h"

struct graph* osh_new_graph(unsigned* offsets, unsigned* adj)
{
  struct graph* g = LOOP_HOST_MALLOC(struct graph, 1);
  g->offsets = offsets;
  g->adj = adj;
  return g;
}

void osh_free_graph(struct graph* g)
{
  if (!g)
    return;
  loop_free(g->offsets);
  loop_free(g->adj);
  loop_host_free(g);
}
