#include "graph.h"
#include "loop.h"

struct graph* new_graph(unsigned* offsets, unsigned* adj)
{
  struct graph* g = loop_malloc(sizeof(*g));
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
  loop_free(g);
}
