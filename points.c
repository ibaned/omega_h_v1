#include "points.h"
#include <stdlib.h>

struct points* new_points(unsigned* offsets, unsigned* adj, double* coords)
{
  struct points* ps = malloc(sizeof(*ps));
  ps->offsets = offsets;
  ps->adj = adj;
  ps->coords = coords;
  return ps;
}

void free_points(struct points* ps)
{
  free(ps->offsets);
  free(ps->adj);
  free(ps->coords);
  free(ps);
}
