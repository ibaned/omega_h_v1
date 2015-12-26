#include "refine_common.h"

#include <assert.h>

#include "arrays.h"
#include "concat_topology.h"
#include "graph.h"
#include "indset.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "quality.h"
#include "refine_class.h"
#include "refine_nodal.h"
#include "refine_qualities.h"
#include "refine_topology.h"
#include "splits_to_domains.h"
#include "tag.h"

static void refine_verts(struct mesh* m, struct mesh* m_out,
    unsigned src_dim, unsigned const* gen_offset_of_srcs)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned nsrcs = mesh_count(m, src_dim);
  unsigned nsplit_srcs = gen_offset_of_srcs[nsrcs];
  unsigned nverts_out = nverts + nsplit_srcs;
  mesh_set_ents(m_out, 0, nverts_out, 0);
  unsigned const* verts_of_srcs = mesh_ask_down(m, src_dim, 0);
  for (unsigned i = 0; i < mesh_count_tags(m, 0); ++i) {
    struct const_tag* t = mesh_get_tag(m, 0, i);
    if (t->type != TAG_F64)
      continue;
    double* gen_vals = refine_nodal(src_dim, nsrcs, verts_of_srcs,
        gen_offset_of_srcs, t->ncomps, t->d.f64);
    double* vals_out = concat_doubles(t->ncomps, t->d.f64, nverts,
        gen_vals, nsplit_srcs);
    loop_free(gen_vals);
    mesh_add_tag(m_out, 0, t->type, t->name, t->ncomps, vals_out);
  }
  refine_class(m, m_out, src_dim, gen_offset_of_srcs);
}

static void refine_ents(struct mesh* m, struct mesh* m_out,
    unsigned src_dim, unsigned const* gen_offset_of_srcs)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nsrcs = mesh_count(m, src_dim);
  unsigned nverts = mesh_count(m, 0);
  unsigned* gen_vert_of_srcs = LOOP_MALLOC(unsigned, nsrcs);
  for (unsigned i = 0; i < nsrcs; ++i)
    if (gen_offset_of_srcs[i] != gen_offset_of_srcs[i + 1])
      gen_vert_of_srcs[i] = nverts + gen_offset_of_srcs[i];
  unsigned* offset_of_same_ents[4] = {0};
  unsigned ngen_ents[4][4] = {{0}};
  unsigned* verts_of_gen_ents[4][4] = {{0}};
  for (unsigned dom_dim = src_dim; dom_dim <= elem_dim; ++dom_dim) {
    if (mesh_get_rep(m) == MESH_REDUCED && dom_dim != elem_dim)
      continue;
    unsigned* offset_of_doms;
    unsigned* direction_of_doms;
    unsigned* vert_of_doms;
    mesh_splits_to_domains(m, dom_dim, src_dim,
        gen_offset_of_srcs, gen_vert_of_srcs,
        &offset_of_doms, &direction_of_doms, &vert_of_doms);
    for (unsigned prod_dim = 1; prod_dim <= dom_dim; ++prod_dim) {
      if (prod_dim != elem_dim) /* reduced testing, fixme */
        continue;
      mesh_refine_topology(m, dom_dim, src_dim, prod_dim,
          offset_of_doms, direction_of_doms, vert_of_doms,
          &ngen_ents[prod_dim][dom_dim], &verts_of_gen_ents[prod_dim][dom_dim]);
    }
    loop_free(direction_of_doms);
    loop_free(vert_of_doms);
    unsigned ndoms = mesh_count(m, dom_dim);
    offset_of_same_ents[dom_dim] = uints_negate_offsets(offset_of_doms, ndoms);
    loop_free(offset_of_doms);
  }
  loop_free(gen_vert_of_srcs);
  for (unsigned ent_dim = 1; ent_dim <= elem_dim; ++ent_dim) {
    if (mesh_get_rep(m) == MESH_REDUCED && ent_dim != elem_dim)
      continue;
    unsigned nents = mesh_count(m, ent_dim);
    unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
    unsigned nents_out;
    unsigned* verts_of_ents_out;
    concat_verts_of_ents(ent_dim,
        nents, verts_of_ents, offset_of_same_ents[ent_dim],
        ngen_ents[ent_dim], verts_of_gen_ents[ent_dim],
        &nents_out, &verts_of_ents_out);
    loop_free(offset_of_same_ents[ent_dim]);
    for (unsigned dom_dim = ent_dim; dom_dim <= elem_dim; ++dom_dim)
      loop_free(verts_of_gen_ents[ent_dim][dom_dim]);
    mesh_set_ents(m_out, ent_dim, nents_out, verts_of_ents_out);
  }
}

unsigned refine_common(
    struct mesh** p_m,
    unsigned src_dim,
    unsigned* candidates,
    double qual_floor,
    unsigned require_better)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  unsigned nsrcs = mesh_count(m, src_dim);
  if (!uints_max(candidates, nsrcs))
    return 0;
  double* src_quals = mesh_refine_qualities(m, src_dim, candidates,
      qual_floor, require_better);
  if (!uints_max(candidates, nsrcs)) {
    loop_free(src_quals);
    return 0;
  }
  unsigned* gen_offset_of_srcs = mesh_indset_offsets(m, src_dim, candidates,
      src_quals);
  loop_free(src_quals);
  struct mesh* m_out = new_mesh(elem_dim);
  mesh_set_rep(m_out, mesh_get_rep(m));
  refine_verts(m, m_out, src_dim, gen_offset_of_srcs);
  refine_ents(m, m_out, src_dim, gen_offset_of_srcs);
  loop_free(gen_offset_of_srcs);
  free_mesh(m);
  *p_m = m_out;
  return 1;
}
