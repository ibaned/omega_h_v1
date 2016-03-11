#include "swap.hpp"

#include <cassert>
#include <cstdio>

#include "arrays.hpp"
#include "comm.hpp"
#include "doubles.hpp"
#include "ghost_mesh.hpp"
#include "indset.hpp"
#include "inherit.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mark.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "parallel_modify.hpp"
#include "quality.hpp"
#include "swap_conserve.hpp"
#include "swap_fit.hpp"
#include "swap_qualities.hpp"
#include "swap_topology.hpp"
#include "tables.hpp"

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
  unsigned ngen_ents = uints_at(gen_offset_of_edges, nedges);
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
  if (mesh_is_parallel(m))
    inherit_globals(m, m_out, ent_dim, same_ent_offsets);
  if (ent_dim == mesh_dim(m)) {
    swap_conserve(m, m_out, gen_offset_of_edges, same_ent_offsets);
    swap_fit(m, m_out, gen_offset_of_edges, same_ent_offsets);
  }
  loop_free(gen_offset_of_edges);
  loop_free(same_ent_offsets);
}

static void swap_interior(
    struct mesh* m)
{
  unsigned const* indset = mesh_find_tag(m, 1, "indset")->d.u32;
  unsigned nedges = mesh_count(m, 1);
  unsigned long total = comm_add_ulong(uints_sum(indset, nedges));
  unsigned const* ring_sizes = mesh_find_tag(m, 1, "ring_size")->d.u32;
  unsigned elem_dim = mesh_dim(m);
  /* vertex handling */
  unsigned nverts = mesh_count(m, 0);
  struct mesh* m_out = new_mesh(elem_dim, mesh_get_rep(m), mesh_is_parallel(m));
  mesh_set_ents(m_out, 0, nverts, 0);
  if (mesh_is_parallel(m))
    mesh_tag_globals(m, 0);
  copy_tags(mesh_tags(m, 0), mesh_tags(m_out, 0), nverts);
  if (mesh_is_parallel(m))
    mesh_parallel_from_tags(m_out, 0);
  /* end vertex handling */
  if (mesh_get_rep(m) == MESH_REDUCED)
    swap_ents(m, m_out, mesh_dim(m), indset, ring_sizes);
  else
    for (unsigned d = 1; d <= mesh_dim(m); ++d)
      swap_ents(m, m_out, d, indset, ring_sizes);
  printf("swapped %10lu %s\n", total, get_ent_name(1, total));
  overwrite_mesh(m, m_out);
}

static unsigned swap_common(
    struct mesh* m,
    unsigned* candidates)
{
  unsigned nedges = mesh_count(m, 1);
  if (!comm_max_uint(uints_max(candidates, nedges))) {
    loop_free(candidates);
    return 0;
  }
  mesh_unmark_boundary(m, 1, candidates);
  if (!comm_max_uint(uints_max(candidates, nedges))) {
    loop_free(candidates);
    return 0;
  }
  double* edge_quals;
  unsigned* ring_sizes;
  mesh_swap_qualities(m, candidates, &edge_quals, &ring_sizes);
  mesh_conform_uints(m, 1, 1, &candidates);
  mesh_conform_doubles(m, 1, 1, &edge_quals);
  mesh_conform_uints(m, 1, 1, &ring_sizes);
  if (!comm_max_uint(uints_max(candidates, nedges))) {
    loop_free(candidates);
    loop_free(edge_quals);
    loop_free(ring_sizes);
    return 0;
  }
  unsigned* indset = mesh_find_indset(m, 1, candidates, edge_quals);
  loop_free(candidates);
  loop_free(edge_quals);
  mesh_add_tag(m, 1, TAG_U32, "indset", 1, indset);
  mesh_add_tag(m, 1, TAG_U32, "ring_size", 1, ring_sizes);
  if (mesh_is_parallel(m)) {
    set_own_ranks_by_indset(m, 1);
    unghost_mesh(m);
  }
  swap_interior(m);
  return 1;
}

unsigned swap_slivers(
    struct mesh* m,
    double good_qual,
    unsigned nlayers)
{
  assert(mesh_dim(m) == 3);
  if (mesh_is_parallel(m)) {
    assert(mesh_get_rep(m) == MESH_FULL);
    mesh_ensure_ghosting(m, 1);
  }
  unsigned elem_dim = mesh_dim(m);
  unsigned* slivers = mesh_mark_slivers(m, good_qual, nlayers);
  unsigned* candidates = mesh_mark_down(m, elem_dim, 1, slivers);
  loop_free(slivers);
  unsigned ret = swap_common(m, candidates);
  return ret;
}
