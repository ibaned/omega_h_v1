#include "refine_common.h"

#include <assert.h>
#include <stdio.h>

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
#include "subset.h"
#include "tables.h"

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
  unsigned* offset_of_doms[4] = {0};
  unsigned* direction_of_doms[4] = {0};
  unsigned* vert_of_doms[4] = {0};
  offset_of_doms[0] = uints_filled(nverts + 1, 0);
  for (unsigned dom_dim = 1; dom_dim <= elem_dim; ++dom_dim) {
    if (!mesh_has_dim(m, dom_dim))
      continue;
    if (dom_dim >= src_dim)
      mesh_splits_to_domains(m, dom_dim, src_dim,
          gen_offset_of_srcs, gen_vert_of_srcs,
          &offset_of_doms[dom_dim], &direction_of_doms[dom_dim],
          &vert_of_doms[dom_dim]);
    else
      offset_of_doms[dom_dim] = uints_filled(mesh_count(m, dom_dim) + 1, 0);
  }
  loop_free(gen_vert_of_srcs);
  unsigned ngen_ents[4][4] = {{0}};
  refined_prod_counts(m, src_dim, offset_of_doms, ngen_ents);
  unsigned gen_offsets[4][5] = {{0}};
  for (unsigned prod_dim = 1; prod_dim <= elem_dim; ++prod_dim) {
    if (mesh_get_rep(m) == MESH_REDUCED && prod_dim != elem_dim)
      continue;
    unsigned nsplit_eqs = offset_of_doms[prod_dim][mesh_count(m, prod_dim)];
    unsigned nsame_eqs = mesh_count(m, prod_dim) - nsplit_eqs;
    gen_offsets[prod_dim][0] = nsame_eqs;
  }
  gen_offsets[0][0] = nverts;
  for (unsigned prod_dim = 0; prod_dim < 4; ++prod_dim) {
    for (unsigned dom_dim = 0; dom_dim < 4; ++dom_dim)
      gen_offsets[prod_dim][dom_dim + 1] =
        gen_offsets[prod_dim][dom_dim] +
        ngen_ents[prod_dim][dom_dim];
  }
  for (unsigned prod_dim = 1; prod_dim <= elem_dim; ++prod_dim) {
    if (mesh_get_rep(m) == MESH_REDUCED && prod_dim != elem_dim)
      continue;
    unsigned nprods = gen_offsets[prod_dim][4];
    unsigned verts_per_prod = the_down_degrees[prod_dim][0];
    unsigned* verts_of_prods = LOOP_MALLOC(unsigned, nprods * verts_per_prod);
    subset_verts_of_doms(m, prod_dim, offset_of_doms[prod_dim], verts_of_prods);
    for (unsigned dom_dim = prod_dim; dom_dim <= elem_dim; ++dom_dim) {
      if (dom_dim < src_dim)
        continue;
      if (mesh_get_rep(m) == MESH_REDUCED && dom_dim != elem_dim)
        continue;
      mesh_refine_topology(m, dom_dim, src_dim, prod_dim,
          offset_of_doms[dom_dim], direction_of_doms[dom_dim],
          vert_of_doms[dom_dim],
          verts_of_prods + 
              (gen_offsets[prod_dim][dom_dim] * verts_per_prod));
    }
    mesh_set_ents(m_out, prod_dim, gen_offsets[prod_dim][4], verts_of_prods);
  }
  refine_class(m, m_out, src_dim, offset_of_doms, ngen_ents, gen_offsets);
  for (unsigned dom_dim = 0; dom_dim <= elem_dim; ++dom_dim) {
    loop_free(offset_of_doms[dom_dim]);
    loop_free(direction_of_doms[dom_dim]);
    loop_free(vert_of_doms[dom_dim]);
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
