#include "reorder.h"

#include "arrays.h"
#include "bfs.h"
#include "ints.h"
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

static unsigned* invert_order(unsigned* o, unsigned n)
{
  unsigned* io = LOOP_HOST_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    io[o[i]] = i;
  return io;
}

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
  unsigned* order = invert_order(queue, nverts);
  loop_host_free(queue);
  return order;
}

LOOP_KERNEL(count_fan_ents,
    unsigned const* vert_num,
    unsigned const* ents_of_verts_offsets,
    unsigned const* ents_of_verts,
    unsigned const* verts_of_ents,
    unsigned verts_per_ent,
    unsigned* nfan_ents)
  unsigned n = vert_num[i];
  unsigned a = ents_of_verts_offsets[i];
  unsigned b = ents_of_verts_offsets[i + 1];
  unsigned nfans = 0;
  for (unsigned j = a; j < b; ++j) {
    unsigned ent = ents_of_verts[j];
    unsigned const* verts_of_ent = verts_of_ents + ent * verts_per_ent;
    unsigned k;
    for (k = 0; k < verts_per_ent; ++k) {
      unsigned oi = verts_of_ent[k];
      unsigned on = vert_num[oi];
      if (on < n)
        break;
    }
    if (k == verts_per_ent)
      ++nfans;
  }
  nfan_ents[i] = nfans;
}

LOOP_KERNEL(fill_fan_ents,
    unsigned const* vert_num,
    unsigned const* ents_of_verts_offsets,
    unsigned const* ents_of_verts,
    unsigned const* verts_of_ents,
    unsigned verts_per_ent,
    unsigned const* fan_ent_offsets,
    unsigned* fan_ents)
  unsigned n = vert_num[i];
  unsigned a = ents_of_verts_offsets[i];
  unsigned b = ents_of_verts_offsets[i + 1];
  unsigned* vert_fans = fan_ents + fan_ent_offsets[i];
  unsigned nfans = 0;
  for (unsigned j = a; j < b; ++j) {
    unsigned ent = ents_of_verts[j];
    unsigned const* verts_of_ent = verts_of_ents + ent * verts_per_ent;
    unsigned k;
    for (k = 0; k < verts_per_ent; ++k) {
      unsigned oi = verts_of_ent[k];
      unsigned on = vert_num[oi];
      if (on < n)
        break;
    }
    if (k == verts_per_ent)
      vert_fans[nfans++] = ent;
  }
}

LOOP_KERNEL(number_fan_ents,
    unsigned const* vert_num,
    unsigned const* fan_ent_offsets,
    unsigned const* fan_ents,
    unsigned const* new_vert_offsets,
    unsigned* ent_num)
  unsigned en = new_vert_offsets[vert_num[i]];
  unsigned a = fan_ent_offsets[i];
  unsigned b = fan_ent_offsets[i + 1];
  for (unsigned j = a; j < b; ++j) {
    unsigned ent = fan_ents[j];
    ent_num[ent] = en++;
  }
}

unsigned* number_ents(struct mesh* m,
    unsigned ent_dim, unsigned const* vert_num)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned nents = mesh_count(m, ent_dim);
  unsigned const* ents_of_verts_offsets = mesh_ask_up(m, 0, ent_dim)->offsets;
  unsigned const* ents_of_verts = mesh_ask_up(m, 0, ent_dim)->adj;
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned* nfan_ents = LOOP_MALLOC(unsigned, nverts);
  LOOP_EXEC(count_fan_ents, nverts, vert_num,
      ents_of_verts_offsets, ents_of_verts,
      verts_of_ents, verts_per_ent, nfan_ents);
  unsigned* fan_ent_offsets = uints_exscan(nfan_ents, nverts);
  assert(uints_at(fan_ent_offsets, nverts) == nents);
  unsigned* fan_ents = LOOP_MALLOC(unsigned, nents);
  LOOP_EXEC(fill_fan_ents, nverts, vert_num,
      ents_of_verts_offsets, ents_of_verts,
      verts_of_ents, verts_per_ent, fan_ent_offsets, fan_ents);
  unsigned* new_vert_nfans = uints_shuffle(nverts, nfan_ents, 1, vert_num);
  loop_free(nfan_ents);
  unsigned* new_vert_offsets = uints_exscan(new_vert_nfans, nverts);
  loop_free(new_vert_nfans);
  assert(uints_at(new_vert_offsets, nverts) == nents);
  unsigned* ent_num = LOOP_MALLOC(unsigned, nents);
  LOOP_EXEC(number_fan_ents, nverts, vert_num, fan_ent_offsets, fan_ents,
      new_vert_offsets, ent_num);
  loop_free(fan_ents);
  loop_free(fan_ent_offsets);
  loop_free(new_vert_offsets);
  return ent_num;
}
