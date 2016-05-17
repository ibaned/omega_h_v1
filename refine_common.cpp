#include "refine_common.hpp"

#include <cassert>
#include <cstdio>

#include "arrays.hpp"
#include "comm.hpp"
#include "ghost_mesh.hpp"
#include "graph.hpp"
#include "indset.hpp"
#include "inherit.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "parallel_modify.hpp"
#include "quality.hpp"
#include "refine_conserve.hpp"
#include "refine_fit.hpp"
#include "refine_nodal.hpp"
#include "refine_qualities.hpp"
#include "refine_topology.hpp"
#include "splits_to_domains.hpp"
#include "subset.hpp"
#include "tables.hpp"

namespace omega_h {

static void refine_verts(struct mesh* m, struct mesh* m_out,
    unsigned src_dim, unsigned const* gen_offset_of_srcs)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned nsrcs = mesh_count(m, src_dim);
  unsigned nsplit_srcs = array_at(gen_offset_of_srcs, nsrcs);
  unsigned nverts_out = nverts + nsplit_srcs;
  mesh_set_ents(m_out, 0, nverts_out, 0);
  unsigned const* verts_of_srcs = mesh_ask_down(m, src_dim, 0);
  for (unsigned i = 0; i < mesh_count_tags(m, 0); ++i) {
    struct const_tag* t = mesh_get_tag(m, 0, i);
    if (t->type != TAG_F64)
      continue;
    double* gen_vals = refine_nodal(src_dim, nsrcs, verts_of_srcs,
        gen_offset_of_srcs, t->ncomps, t->d.f64);
    double* vals_out = concat_arrays(t->d.f64, gen_vals,
        nverts, nsplit_srcs, t->ncomps);
    loop_free(gen_vals);
    mesh_add_tag(m_out, 0, t->type, t->name, t->ncomps, vals_out);
  }
  if (mesh_is_parallel(m)) {
    unsigned* offsets = uints_linear(nverts + 1, 1);
    inherit_globals(m, m_out, 0, offsets);
    loop_free(offsets);
  }
}

LOOP_KERNEL(gen_vert_of_src_kern,
    unsigned const* gen_offset_of_srcs,
    unsigned nverts,
    unsigned* gen_vert_of_srcs)
  if (gen_offset_of_srcs[i] != gen_offset_of_srcs[i + 1])
    gen_vert_of_srcs[i] = nverts + gen_offset_of_srcs[i];
}

static void refine_ents(struct mesh* m, struct mesh* m_out,
    unsigned src_dim, unsigned const* gen_offset_of_srcs)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nsrcs = mesh_count(m, src_dim);
  unsigned nverts = mesh_count(m, 0);
  unsigned* gen_vert_of_srcs = LOOP_MALLOC(unsigned, nsrcs);
  LOOP_EXEC(gen_vert_of_src_kern, nsrcs,
      gen_offset_of_srcs, nverts, gen_vert_of_srcs);
  unsigned* offset_of_doms[4] = {0};
  unsigned* direction_of_doms[4] = {0};
  unsigned* vert_of_doms[4] = {0};
  offset_of_doms[0] = filled_array<unsigned>(nverts + 1, 0);
  for (unsigned dom_dim = 1; dom_dim <= elem_dim; ++dom_dim) {
    if (!mesh_has_dim(m, dom_dim))
      continue;
    if (dom_dim >= src_dim)
      mesh_splits_to_domains(m, dom_dim, src_dim,
          gen_offset_of_srcs, gen_vert_of_srcs,
          &offset_of_doms[dom_dim], &direction_of_doms[dom_dim],
          &vert_of_doms[dom_dim]);
    else
      offset_of_doms[dom_dim] = filled_array<unsigned>(
          mesh_count(m, dom_dim) + 1, 0);
  }
  loop_free(gen_vert_of_srcs);
  for (unsigned prod_dim = 0; prod_dim <= elem_dim; ++prod_dim) {
    if (mesh_get_rep(m) == MESH_REDUCED &&
        (0 < prod_dim && prod_dim < elem_dim))
      continue;
    unsigned ndoms[4];
    unsigned* prods_of_doms_offsets[4];
    unsigned ngen_offsets[5];
    setup_refine(m, src_dim, prod_dim, offset_of_doms,
        ndoms, prods_of_doms_offsets, ngen_offsets);
    if (prod_dim) {
      unsigned nprods_out = ngen_offsets[4];
      unsigned verts_per_prod = the_down_degrees[prod_dim][0];
      unsigned* verts_of_prods_out = LOOP_MALLOC(
          unsigned, nprods_out * verts_per_prod);
      expand_into(verts_of_prods_out, mesh_ask_down(m, prod_dim, 0),
          prods_of_doms_offsets[0], ndoms[0], verts_per_prod);
      for (unsigned dom_dim = 1; dom_dim <= elem_dim; ++dom_dim)
        if (ndoms[dom_dim] && prods_of_doms_offsets[dom_dim])
          mesh_refine_topology(m, dom_dim, src_dim, prod_dim,
              offset_of_doms[dom_dim], direction_of_doms[dom_dim],
              vert_of_doms[dom_dim],
              verts_of_prods_out + (ngen_offsets[dom_dim] * verts_per_prod));
      mesh_set_ents(m_out, prod_dim, nprods_out, verts_of_prods_out);
      if (mesh_is_parallel(m))
        inherit_globals(m, m_out, prod_dim, prods_of_doms_offsets[0]);
    }
    inherit_class(m, m_out, prod_dim, ndoms, prods_of_doms_offsets);
    if (prod_dim == elem_dim) {
      refine_conserve(m, m_out, ndoms, prods_of_doms_offsets);
      refine_fit(m, m_out, ndoms, prods_of_doms_offsets);
    }
    for (unsigned i = 0; i < 4; ++i)
      loop_free(prods_of_doms_offsets[i]);
  }
  for (unsigned dom_dim = 0; dom_dim <= elem_dim; ++dom_dim) {
    loop_free(offset_of_doms[dom_dim]);
    loop_free(direction_of_doms[dom_dim]);
    loop_free(vert_of_doms[dom_dim]);
  }
}

static unsigned choose_refinement_indset(
    struct mesh* m,
    unsigned src_dim,
    double qual_floor,
    unsigned require_better)
{
  if (mesh_is_parallel(m))
    mesh_ensure_ghosting(m, 1);
  unsigned nsrcs = mesh_count(m, src_dim);
  unsigned const* candidates = mesh_find_tag(m, src_dim, "candidate")->d.u32;
  if (!comm_max_uint(uints_max(candidates, nsrcs))) {
    mesh_free_tag(m, src_dim, "candidate");
    fprintf(stderr, "no candidates at all\n");
    return 0;
  }
  unsigned* good_candidates = copy_array(candidates, nsrcs);
  mesh_free_tag(m, src_dim, "candidate");
  double* src_quals = mesh_refine_qualities(m, src_dim, &good_candidates,
      qual_floor, require_better);
  if (!comm_max_uint(uints_max(good_candidates, nsrcs))) {
    loop_free(good_candidates);
    loop_free(src_quals);
    fprintf(stderr, "quality checks failed\n");
    return 0;
  }
  unsigned* indset = mesh_find_indset(m, src_dim,
      good_candidates, src_quals);
  loop_free(good_candidates);
  loop_free(src_quals);
  mesh_add_tag(m, src_dim, TAG_U32, "indset", 1, indset);
  return 1;
}

unsigned refine_common(
    struct mesh* m,
    unsigned src_dim,
    double qual_floor,
    unsigned require_better)
{
  if (mesh_is_parallel(m))
    assert(mesh_get_rep(m) == MESH_FULL);
  if (!choose_refinement_indset(m, src_dim, qual_floor, require_better)) {
    fprintf(stderr, "refine_common returning zero due to indset\n");
    return 0;
  }
  if (mesh_is_parallel(m)) {
    set_own_ranks_by_indset(m, src_dim);
    unghost_mesh(m);
  }
  unsigned const* indset = mesh_find_tag(m, src_dim, "indset")->d.u32;
  unsigned nsrcs = mesh_count(m, src_dim);
  unsigned long total = comm_add_ulong(uints_sum(indset, nsrcs));
  unsigned* gen_offset_of_srcs = uints_exscan(indset, nsrcs);
  mesh_free_tag(m, src_dim, "indset");
  struct mesh* m_out = new_mesh(mesh_dim(m),
      mesh_get_rep(m), mesh_is_parallel(m));
  refine_verts(m, m_out, src_dim, gen_offset_of_srcs);
  refine_ents(m, m_out, src_dim, gen_offset_of_srcs);
  loop_free(gen_offset_of_srcs);
  if (comm_rank() == 0)
    printf("split %10lu %s\n", total, get_ent_name(src_dim, total));
  overwrite_mesh(m, m_out);
  fprintf(stderr, "refine_common returning one after overwrite\n");
  return 1;
}

}
