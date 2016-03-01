#include "reorder.h"

#include "arrays.h"
#include "bfs.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "tables.h"

static void get_host_graph(struct mesh* m,
    unsigned** p_offsets,
    unsigned** p_adj)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned* offsets = uints_to_host(
      mesh_ask_star(m, 0, 1)->offsets, nverts + 1);
  unsigned nadj = offsets[nverts];
  unsigned* adj = uints_to_host(
      mesh_ask_star(m, 0, 1)->adj, nadj);
  *p_offsets = offsets;
  *p_adj = adj;
}

static unsigned* compute_boundary_depth(struct mesh* m,
    unsigned const* offsets, unsigned const* adj)
{
  unsigned dim = mesh_dim(m);
  unsigned nverts = mesh_count(m, 0);
  unsigned* bdry_sides = mesh_mark_part_boundary(m);
  unsigned* bdry_verts = mesh_mark_down(m, dim - 1, 0, bdry_sides);
  loop_free(bdry_sides);
  unsigned* host_bdry_verts = uints_to_host(bdry_verts, nverts);
  loop_free(bdry_verts);
  unsigned* queue = LOOP_HOST_MALLOC(unsigned, nverts);
  unsigned* depth = LOOP_HOST_MALLOC(unsigned, nverts);
  unsigned begin = 0;
  unsigned end = 0;
  for (unsigned i = 0; i < nverts; ++i)
    if (host_bdry_verts[i]) {
      queue[end++] = i;
      depth[i] = 0;
    } else {
      depth[i] = INVALID;
    }
  loop_host_free(host_bdry_verts);
  bfs_continue(queue, &begin, &end, offsets, adj, 0, depth);
  loop_host_free(queue);
  return depth;
}

static void enqueue_component_maxima(
    unsigned nverts,
    unsigned const* depth, unsigned const* comp,
    unsigned* queue, unsigned* p_end)
{
  for (unsigned c = 0; c < nverts; ++c) {
    unsigned max_v = INVALID;
    unsigned max_depth = INVALID;
    unsigned exists = 0;
    for (unsigned v = 0; v < nverts; ++v) {
      if (comp[v] != c)
        continue;
      exists = 1;
      if (max_v == INVALID || depth[v] > max_depth) {
        max_v = v;
        max_depth = depth[v];
      }
    }
    if (!exists)
      break;
    assert(max_v != INVALID);
    queue[(*p_end)++] = max_v;
  }
}

/* returns a map from NEW indices to OLD indices */
unsigned* compute_ordering(struct mesh* m)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned* offsets;
  unsigned* adj;
  get_host_graph(m, &offsets, &adj);
  unsigned* depth = compute_boundary_depth(m, offsets, adj);
  unsigned* comp = LOOP_HOST_MALLOC(unsigned, nverts);
  connected_components(nverts, offsets, adj, comp);
  unsigned* queue = LOOP_HOST_MALLOC(unsigned, nverts);
  unsigned end = 0;
  enqueue_component_maxima(nverts, depth, comp, queue, &end);
  unsigned ncomps = end;
  unsigned begin = 0;
  unsigned* radius = depth;
  for (unsigned i = 0; i < nverts; ++i)
    radius[i] = INVALID;
  for (unsigned i = 0; i < end; ++i)
    radius[queue[i]] = 0;
  bfs_continue(queue, &begin, &end, offsets, adj, 0, radius);
  assert(end == nverts);
  /* we could do all the components from their boundary seeds
     at once and then sort by component first and layer second.
     we'll avoid the need to bring a sort into this by doing
     one component at a time. */
  unsigned* comp_seeds = LOOP_HOST_MALLOC(unsigned, ncomps);
  end = 0;
  enqueue_component_maxima(nverts, radius, comp, comp_seeds, &end);
  assert(end == ncomps);
  loop_host_free(comp);
  unsigned* layer = radius;
  for (unsigned i = 0; i < nverts; ++i)
    layer[i] = INVALID;
  begin = 0;
  end = 0;
  for (unsigned c = 0; c < ncomps; ++c) {
    queue[end++] = comp_seeds[c];
    layer[comp_seeds[c]] = 0;
    bfs_continue(queue, &begin, &end, offsets, adj, 0, layer);
  }
  assert(end == nverts);
  loop_host_free(comp_seeds);
  loop_host_free(offsets);
  loop_host_free(adj);
  loop_host_free(layer);
  return queue;
}
