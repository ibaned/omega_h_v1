#include "up_from_down.h"
#include "tables.h"
#include "ints.h"
#include <stdlib.h>

struct up_from_down get_up_from_down(
    unsigned up_dim,
    unsigned down_dim,
    unsigned nup,
    unsigned ndown,
    unsigned* edges)
{
  struct up_from_down out;
  unsigned down_degree = the_down_degrees[up_dim][down_dim];
  unsigned nedges = nup * down_degree;
  unsigned* up_degrees = malloc(sizeof(unsigned) * ndown);
  ints_zero(up_degrees, ndown);
  for (unsigned i = 0; i < nedges; ++i)
    up_degrees[edges[i]]++;
  out.offsets = ints_exscan(up_degrees, ndown);
  out.edges = malloc(sizeof(unsigned) * nedges);
  out.directions = malloc(sizeof(unsigned) * nedges);
  ints_zero(up_degrees, ndown);
  for (unsigned i = 0; i < nedges; ++i) {
    unsigned up = i / down_degree;
    unsigned direction = i % down_degree;
    unsigned down = edges[i];
    unsigned o = out.offsets[down];
    unsigned j = up_degrees[down]++;
    out.edges[o + j] = up;
    out.directions[o + j] = direction;
  }
  free(up_degrees);
  return out;
}
