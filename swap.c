#include "swap.h"

#include <assert.h>

#include "comm.h"
#include "copy_tags.h"
#include "doubles.h"
#include "ghost_mesh.h"
#include "indset.h"
#include "inherit.h"
#include "ints.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "parallel_modify.h"
#include "quality.h"
#include "swap_conserve.h"
#include "swap_qualities.h"
#include "swap_topology.h"
#include "tables.h"

static void swap_ents(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ent_dim,
    unsigned const* indset,
    unsigned const* ring_sizes)
{
  unsigned nedges = mesh_count(m, 1);
  unsigned* gen_offset_of_edges = get_swap_topology_offsets(
      ent_dim, nedges, indset, ring_sizes);
  unsigned ngen_ents = gen_offset_of_edges[nedges];
  unsigned* verts_of_gen_ents = mesh_swap_topology(m, ent_dim,
      indset, gen_offset_of_edges);
  unsigned* old_ents = mesh_mark_up(m, 1, ent_dim, indset);
  unsigned nents = mesh_count(m, ent_dim);
  unsigned* same_ents = uints_negate(old_ents, nents);
  loop_free(old_ents);
  unsigned* same_ent_offsets = uints_exscan(same_ents, nents);
  loop_free(same_ents);
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned nents_out;
  unsigned* verts_of_ents_out;
  concat_verts_of_ents(ent_dim, nents, ngen_ents, verts_of_ents,
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
  if (ent_dim == mesh_dim(m))
    swap_conserve(m, m_out, gen_offset_of_edges, same_ent_offsets);
  loop_free(gen_offset_of_edges);
  loop_free(same_ent_offsets);
}

static void swap_interior(
    struct mesh** p_m)
{
  struct mesh* m = *p_m;
  unsigned const* indset = mesh_find_tag(m, 1, "indset")->d.u32;
  unsigned const* ring_sizes = mesh_find_tag(m, 1, "ring_size")->d.u32;
  unsigned elem_dim = mesh_dim(m);
  /* vertex handling */
  unsigned nverts = mesh_count(m, 0);
  struct mesh* m_out = new_mesh(elem_dim, mesh_get_rep(m), mesh_is_parallel(m));
  mesh_set_ents(m_out, 0, nverts, 0);
  copy_tags(mesh_tags(m, 0), mesh_tags(m_out, 0), nverts);
  /* end vertex handling */
  if (mesh_get_rep(m) == MESH_REDUCED)
    swap_ents(m, m_out, mesh_dim(m), indset, ring_sizes);
  else
    for (unsigned d = 1; d <= mesh_dim(m); ++d)
      swap_ents(m, m_out, d, indset, ring_sizes);
  free_mesh(m);
  *p_m = m_out;
}

static unsigned swap_common(
    struct mesh** p_m,
    unsigned* candidates)
{
  struct mesh* m = *p_m;
  unsigned nedges = mesh_count(m, 1);
  if (!comm_max_uint(uints_max(candidates, nedges)))
    return 0;
  mesh_unmark_boundary(m, 1, candidates);
  if (!comm_max_uint(uints_max(candidates, nedges)))
    return 0;
  double* edge_quals;
  unsigned* ring_sizes;
  mesh_swap_qualities(m, candidates, &edge_quals, &ring_sizes);
  if (!comm_max_uint(uints_max(candidates, nedges))) {
    loop_free(edge_quals);
    loop_free(ring_sizes);
    return 0;
  }
  unsigned* indset = mesh_find_indset(m, 1, candidates, edge_quals);
  loop_free(edge_quals);
  mesh_add_tag(m, 1, TAG_U32, "indset", 1, indset);
  mesh_add_tag(m, 1, TAG_U32, "ring_size", 1, ring_sizes);
  if (mesh_is_parallel(*p_m)) {
    set_own_ranks_by_indset(*p_m, 1);
    unghost_mesh(p_m);
  }
  swap_interior(p_m);
  return 1;
}

unsigned swap_slivers(
    struct mesh** p_m,
    double good_qual,
    unsigned nlayers)
{
  assert(mesh_dim(*p_m) == 3);
  if (mesh_is_parallel(*p_m)) {
    assert(mesh_get_rep(*p_m) == MESH_FULL);
    mesh_ensure_ghosting(p_m, 1);
  }
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned* slivers = mesh_mark_slivers(m, good_qual, nlayers);
  unsigned* candidates = mesh_mark_down(m, elem_dim, 1, slivers);
  loop_free(slivers);
  unsigned ret = swap_common(p_m, candidates);
  loop_free(candidates);
  return ret;
}
