#include "subset.h"
#include "ints.h"
#include <stdlib.h>

unsigned* ints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* set_offsets)
{
  unsigned nsub = set_offsets[n];
  unsigned* out = malloc(sizeof(unsigned) * nsub * width);
  for (unsigned i = 0; i < n; ++i) {
    if (set_offsets[i] == set_offsets[i + 1])
      continue;
    for (unsigned j = 0; j < width; ++j)
      out[set_offsets[i] * width + j] = a[i * width + j];
  }
  return out;
}

struct subgraph induce_subgraph(
    unsigned nverts,
    unsigned const* offsets,
    unsigned const* edges,
    unsigned const* set_offsets)
{
  unsigned* degrees = malloc(sizeof(unsigned) * nverts);
  for (unsigned i = 0; i < nverts; ++i)
    degrees[i] = offsets[i + 1] - offsets[i];
  unsigned* sub_degrees = ints_subset(nverts, 1, degrees, set_offsets);
  free(degrees);
  struct subgraph out;
  out.nverts = set_offsets[nverts];
  out.offsets = ints_exscan(sub_degrees, out.nverts);
  free(sub_degrees);
  unsigned nsub_edges = out.offsets[out.nverts];
  out.edges = malloc(sizeof(unsigned) * nsub_edges);
  for (unsigned i = 0; i < nverts; ++i) {
    if (set_offsets[i] == set_offsets[i + 1])
      continue;
    unsigned sub_vert = set_offsets[i];
    unsigned first_edge = offsets[i];
    unsigned end_edge = offsets[i + 1];
    unsigned first_sub_edge = out.offsets[sub_vert];
    for (unsigned j = first_edge; j < end_edge; ++j)
      out.edges[first_sub_edge + (j - first_edge)] = edges[j];
  }
  return out;
}
