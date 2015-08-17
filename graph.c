#include "graph.h"
#include <stdlib.h>

struct graph* new_graph(unsigned* offsets, unsigned* adj)
{
  struct graph* g = malloc(sizeof(*g));
  g->offsets = offsets;
  g->adj = adj;
  return g;
}

void free_graph(struct graph* g)
{
  if (!g)
    return;
  free(g->offsets);
  free(g->adj);
  free(g);
}
