#include "swap_common.h"

#include <assert.h>

#include "copy_tags.h"
#include "doubles.h"
#include "graph.h"
#include "indset.h"
#include "inherit.h"
#include "ints.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "quality.h"
#include "swap_qualities.h"
#include "swap_topology.h"
#include "tables.h"

static void swap_ents(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ent_dim,
    unsigned const* indset,
    unsigned const* ring_sizes,
    unsigned const* edge_codes)
{
  unsigned nedges = mesh_count(m, 1);
  unsigned* gen_offset_of_edges = get_swap_topology_offsets(
      ent_dim, nedges, indset, ring_sizes);
  unsigned ngen_ents = gen_offset_of_edges[nedges];
  unsigned* verts_of_gen_ents = mesh_swap_topology(m, ent_dim,
      indset, gen_offset_of_edges, edge_codes);
  unsigned* old_ents = mesh_mark_up(m, 1, ent_dim, indset);
  unsigned nents = mesh_count(m, ent_dim);
  unsigned* same_ents = uints_negate(old_ents, nents);
  loop_free(old_ents);
  unsigned* same_ent_offsets = uints_exscan(same_ents, nents);
  loop_free(same_ents);
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned nents_out;
  unsigned* verts_of_ents_out;
  concat_verts_of_elems(ent_dim, nents, ngen_ents, verts_of_ents,
      same_ent_offsets, verts_of_gen_ents,
      &nents_out, &verts_of_ents_out);
  loop_free(verts_of_gen_ents);
  mesh_set_ents(m_out, ent_dim, nents_out, verts_of_ents_out);
  if (mesh_get_rep(m) == MESH_FULL) {
    unsigned ndoms[4];
    unsigned* prods_of_doms_offsets[4];
    setup_swap(m, ent_dim, gen_offset_of_edges, same_ent_offsets,
        ndoms, prods_of_doms_offsets);
    inherit_class(m, m_out, ent_dim, ndoms, prods_of_doms_offsets);
  }
  loop_free(gen_offset_of_edges);
  loop_free(same_ent_offsets);
}

unsigned swap_common(
    struct mesh** p_m,
    unsigned* candidates)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  assert(elem_dim == 3);
  unsigned nedges = mesh_count(m, 1);
  if (!uints_max(candidates, nedges))
    return 0;
  mesh_unmark_boundary(m, 1, candidates);
  if (!uints_max(candidates, nedges))
    return 0;
  double* edge_quals;
  unsigned* edge_codes;
  unsigned* ring_sizes;
  mesh_swap_qualities(m, candidates,
      &edge_quals, &edge_codes, &ring_sizes);
  if (!uints_max(candidates, nedges)) {
    loop_free(edge_quals);
    loop_free(edge_codes);
    loop_free(ring_sizes);
    return 0;
  }
  /* vertex handling */
  unsigned nverts = mesh_count(m, 0);
  struct mesh* m_out = new_mesh(elem_dim, mesh_get_rep(m), 0);
  mesh_set_ents(m_out, 0, nverts, 0);
  copy_tags(mesh_tags(m, 0), mesh_tags(m_out, 0), nverts);
  /* end vertex handling */
  unsigned* indset = mesh_find_indset(m, 1, candidates, edge_quals);
  loop_free(edge_quals);
  if (mesh_get_rep(m) == MESH_REDUCED)
    swap_ents(m, m_out, mesh_dim(m), indset, ring_sizes, edge_codes);
  else
    for (unsigned d = 1; d <= mesh_dim(m); ++d)
      swap_ents(m, m_out, d, indset, ring_sizes, edge_codes);
  loop_free(indset);
  loop_free(edge_codes);
  loop_free(ring_sizes);
  free_mesh(m);
  *p_m = m_out;
  return 1;
}
