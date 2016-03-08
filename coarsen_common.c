#include "coarsen_common.h"

#include <stdio.h>

#include "arrays.h"
#include "check_collapse_class.h"
#include "coarsen_conserve.h"
#include "coarsen_fit.h"
#include "coarsen_qualities.h"
#include "coarsen_topology.h"
#include "collapse_codes.h"
#include "collapses_to_ents.h"
#include "collapses_to_verts.h"
#include "comm.h"
#include "ghost_mesh.h"
#include "indset.h"
#include "inherit.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "parallel_modify.h"
#include "quality.h"
#include "subset.h"
#include "tables.h"

LOOP_KERNEL(remap_conn,
    unsigned const* offset_of_same_verts,
    unsigned* verts_of_ents_out)
  verts_of_ents_out[i] = offset_of_same_verts[verts_of_ents_out[i]];
}

static void coarsen_ents(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ent_dim,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned const* offset_of_same_verts,
    unsigned const* fused_ents,
    unsigned** p_fused_sides)
{
  unsigned nents = mesh_count(m, ent_dim);
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned* gen_offset_of_ents;
  unsigned* gen_vert_of_ents;
  unsigned* gen_direction_of_ents;
  unsigned* offset_of_same_ents;
  collapses_to_ents(m, ent_dim,
      gen_offset_of_verts, gen_vert_of_verts, fused_ents,
      &gen_offset_of_ents, &gen_vert_of_ents,
      &gen_direction_of_ents, &offset_of_same_ents, p_fused_sides);
  unsigned ngen_ents;
  unsigned* verts_of_gen_ents;
  coarsen_topology(ent_dim, nents, verts_of_ents, gen_offset_of_ents,
      gen_vert_of_ents, gen_direction_of_ents, &ngen_ents,
      &verts_of_gen_ents);
  loop_free(gen_direction_of_ents);
  loop_free(gen_vert_of_ents);
  unsigned nents_out;
  unsigned* verts_of_ents_out;
  concat_verts_of_ents(ent_dim, nents, ngen_ents, verts_of_ents,
      offset_of_same_ents, verts_of_gen_ents,
      &nents_out, &verts_of_ents_out);
  loop_free(verts_of_gen_ents);
  /* remap new connectivity to account for vertex removal */
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  LOOP_EXEC(remap_conn, nents_out * verts_per_ent,
      offset_of_same_verts,
      verts_of_ents_out);
  mesh_set_ents(m_out, ent_dim, nents_out, verts_of_ents_out);
  unsigned ndoms[4];
  unsigned* prods_of_doms_offsets[4];
  setup_coarsen(m, ent_dim, gen_offset_of_ents, offset_of_same_ents,
      ndoms, prods_of_doms_offsets);
  inherit_class(m, m_out, ent_dim, ndoms, prods_of_doms_offsets);
  if (mesh_is_parallel(m))
    inherit_globals(m, m_out, ent_dim, offset_of_same_ents);
  if (ent_dim == mesh_dim(m)) {
    coarsen_conserve(m, m_out, gen_offset_of_verts, gen_offset_of_ents,
        offset_of_same_ents);
    coarsen_fit(m, m_out, gen_offset_of_verts, gen_offset_of_ents,
        offset_of_same_ents);
  }
  loop_free(gen_offset_of_ents);
  loop_free(offset_of_same_ents);
}

static void coarsen_all_ents(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned const* offset_of_same_verts)
{
  unsigned* fused_sides[4] = {0};
  for (unsigned dd = 0; dd < mesh_dim(m); ++dd) {
    unsigned d = mesh_dim(m) - dd;
    coarsen_ents(m, m_out, d, gen_offset_of_verts,
        gen_vert_of_verts, offset_of_same_verts,
        fused_sides[d], &fused_sides[d - 1]);
  }
  for (unsigned d = 0; d <= mesh_dim(m); ++d)
    loop_free(fused_sides[d]);
}

static unsigned check_coarsen_noop(struct mesh* m)
{
  unsigned nedges = mesh_count(m, 1);
  unsigned const* col_codes = mesh_find_tag(m, 1, "col_codes")->d.u32;
  if (comm_max_uint(uints_max(col_codes, nedges)) == DONT_COLLAPSE) {
    mesh_free_tag(m, 1, "col_codes");
    return 0;
  }
  return 1;
}

static unsigned check_coarsen_class(struct mesh* m)
{
  /* right now this assumes we're doing the simple
     check_collapse_class that only looks at the edge
     closure classification. if it gets more advanced,
     add ghosting and synchronization to this function */
  unsigned const* col_codes_in = mesh_find_tag(m, 1, "col_codes")->d.u32;
  unsigned nedges = mesh_count(m, 1);
  unsigned* col_codes = uints_copy(col_codes_in, nedges);
  mesh_free_tag(m, 1, "col_codes");
  check_collapse_class(m, col_codes);
  if (comm_max_uint(uints_max(col_codes, nedges)) == DONT_COLLAPSE) {
    loop_free(col_codes);
    return 0;
  }
  mesh_add_tag(m, 1, TAG_U32, "col_codes", 1, col_codes);
  return 1;
}

static unsigned check_coarsen_quality(
    struct mesh* m,
    double quality_floor,
    unsigned require_better)
{
  if (mesh_is_parallel(m))
    mesh_ensure_ghosting(m, 1);
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned nedges = mesh_count(m, 1);
  unsigned const* col_codes_in = mesh_find_tag(m, 1, "col_codes")->d.u32;
  unsigned* col_codes = uints_copy(col_codes_in, nedges);
  mesh_free_tag(m, 1, "col_codes");
  double const* coords = mesh_find_tag(m, 0, "coordinates")->d.f64;
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* elems_of_verts_offsets =
    mesh_ask_up(m, 0, elem_dim)->offsets;
  unsigned const* elems_of_verts =
    mesh_ask_up(m, 0, elem_dim)->adj;
  unsigned const* elems_of_verts_directions =
    mesh_ask_up(m, 0, elem_dim)->directions;
  double* elem_quals = 0;
  if (require_better)
    elem_quals = element_qualities(elem_dim, nelems, verts_of_elems, coords);
  double* quals_of_edges = coarsen_qualities(elem_dim, nedges, col_codes,
      verts_of_elems, verts_of_edges,
      elems_of_verts_offsets, elems_of_verts, elems_of_verts_directions,
      coords, quality_floor, elem_quals);
  loop_free(elem_quals);
  mesh_conform_uints(m, 1, 1, &col_codes);
  if (comm_max_uint(uints_max(col_codes, nedges)) == DONT_COLLAPSE) {
    loop_free(col_codes);
    loop_free(quals_of_edges);
    return 0;
  }
  mesh_conform_doubles(m, 1, 2, &quals_of_edges);
  mesh_add_tag(m, 1, TAG_U32, "col_codes", 1, col_codes);
  mesh_add_tag(m, 1, TAG_F64, "col_quals", 2, quals_of_edges);
  return 1;
}

static void setup_coarsen_indset(struct mesh* m)
{
  valid_collapses_to_verts(m);
  unsigned const* candidates = mesh_find_tag(m, 0, "candidates")->d.u32;
  double const* quals = mesh_find_tag(m, 0, "col_qual")->d.f64;
  unsigned* indset = mesh_find_indset(m, 0, candidates, quals);
  mesh_free_tag(m, 0, "candidates");
  mesh_add_tag(m, 0, TAG_U32, "indset", 1, indset);
}

static void coarsen_interior(struct mesh* m)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nverts = mesh_count(m, 0);
  unsigned* gen_vert_of_verts = collapsing_vertex_destinations(m);
  mesh_free_tag(m, 1, "col_codes");
  mesh_free_tag(m, 1, "col_quals");
  mesh_free_tag(m, 0, "col_qual");
  unsigned const* indset = mesh_find_tag(m, 0, "indset")->d.u32;
  unsigned long total = comm_add_ulong(uints_sum(indset, nverts));
  unsigned* gen_offset_of_verts = uints_exscan(indset, nverts);
  mesh_free_tag(m, 0, "indset");
  unsigned* offset_of_same_verts = uints_negate_offsets(
      gen_offset_of_verts, nverts);
  unsigned nverts_out = uints_at(offset_of_same_verts, nverts);
  struct mesh* m_out = new_mesh(elem_dim, mesh_get_rep(m), mesh_is_parallel(m));
  mesh_set_ents(m_out, 0, nverts_out, 0);
  tags_subset(m, m_out, 0, offset_of_same_verts);
  if (mesh_is_parallel(m))
    inherit_globals(m, m_out, 0, offset_of_same_verts);
  coarsen_all_ents(m, m_out, gen_offset_of_verts, gen_vert_of_verts,
      offset_of_same_verts);
  loop_free(gen_vert_of_verts);
  loop_free(gen_offset_of_verts);
  loop_free(offset_of_same_verts);
  if (comm_rank() == 0)
    printf("collapsed %10lu %s\n", total, get_ent_name(1, total));
  overwrite_mesh(m, m_out);
}

unsigned coarsen_common(
    struct mesh* m,
    double quality_floor,
    unsigned require_better)
{
  if (!check_coarsen_noop(m))
    return 0;
  if (!check_coarsen_class(m))
    return 0;
  if (!check_coarsen_quality(m, quality_floor, require_better))
    return 0;
  setup_coarsen_indset(m);
  if (mesh_is_parallel(m)) {
    set_own_ranks_by_indset(m, 0);
    unghost_mesh(m);
  }
  coarsen_interior(m);
  return 1;
}
