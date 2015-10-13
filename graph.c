#include "graph.h"

#include "loop.h"

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
  LOOP_FREE(g->offsets);
  LOOP_FREE(g->adj);
  LOOP_HOST_FREE(g);
}
