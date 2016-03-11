#include "bfs.hpp"

#include <cassert>

#include "loop.hpp"
#include "tables.hpp"

/* none of the functions in this file
   are parallelized via the LOOP mechanism.
   this means during accelerator runs they
   will have to be executed on the host.
   the reasoning is that BFS is somewhat
   inherently serial (#layers), and parallelizing
   each layer seems like more trouble than its worth */

void bfs_continue(
    unsigned* queue,
    unsigned* p_begin,
    unsigned* p_end,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned* comp,
    unsigned* layer)
{
  unsigned begin = *p_begin;
  unsigned end = *p_end;
  while (begin != end) {
    unsigned u = queue[begin++];
    unsigned a = offsets[u];
    unsigned b = offsets[u + 1];
    for (unsigned j = a; j < b; ++j) {
      unsigned v = adj[j];
      if (layer[v] == INVALID) {
        queue[end++] = v;
        layer[v] = layer[u] + 1;
        if (comp)
          comp[v] = comp[u];
      }
    }
  }
  *p_begin = begin;
  *p_end = end;
}

void bfs_from(
    unsigned start,
    unsigned the_comp,
    unsigned* queue,
    unsigned* p_begin,
    unsigned* p_end,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned* comp,
    unsigned* layer)
{
  unsigned begin = *p_begin;
  unsigned end = *p_end;
  queue[end++] = start;
  layer[start] = 0;
  if (comp)
    comp[start] = the_comp;
  bfs_continue(queue, &begin, &end,
      offsets, adj, comp, layer);
  *p_begin = begin;
  *p_end = end;
}

void bfs_full(
    unsigned n,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned* comp,
    unsigned* layer,
    unsigned* sorted)
{
  for (unsigned i = 0; i < n; ++i)
    layer[i] = INVALID;
  unsigned begin = 0;
  unsigned end = 0;
  unsigned the_comp = 0;
  for (unsigned i = 0; i < n; ++i)
    if (layer[i] == INVALID)
      bfs_from(i, the_comp++, sorted, &begin, &end,
          offsets, adj, comp, layer);
  assert(end == n);
}

void connected_components(
    unsigned n,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned* comp)
{
  unsigned* layer = LOOP_HOST_MALLOC(unsigned, n);
  unsigned* sorted = LOOP_HOST_MALLOC(unsigned, n);
  bfs_full(n, offsets, adj, comp, layer, sorted);
  loop_host_free(layer);
  loop_host_free(sorted);
}

