#include "refine_common.h"

#include <assert.h>

#include "arrays.h"
#include "graph.h"
#include "indset.h"
#include "inherit.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "quality.h"
#include "refine_conserve.h"
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
      unsigned* verts_of_prods_out = LOOP_MALLOC(unsigned, nprods_out * verts_per_prod);
      uints_expand_into(ndoms[0], verts_per_prod, mesh_ask_down(m, prod_dim, 0),
          prods_of_doms_offsets[0], verts_of_prods_out);
      for (unsigned dom_dim = 1; dom_dim <= elem_dim; ++dom_dim)
        if (ndoms[dom_dim] && prods_of_doms_offsets[dom_dim])
          mesh_refine_topology(m, dom_dim, src_dim, prod_dim,
              offset_of_doms[dom_dim], direction_of_doms[dom_dim],
              vert_of_doms[dom_dim],
              verts_of_prods_out + (ngen_offsets[dom_dim] * verts_per_prod));
      mesh_set_ents(m_out, prod_dim, nprods_out, verts_of_prods_out);
    }
    inherit_class(m, m_out, prod_dim, ndoms, prods_of_doms_offsets);
    if (prod_dim == elem_dim)
      refine_conserve(m, m_out, ndoms, prods_of_doms_offsets);
    for (unsigned i = 0; i < 4; ++i)
      loop_free(prods_of_doms_offsets[i]);
  }
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
  if (mesh_is_parallel(m))
    assert(mesh_ghost_layers(m) == 1);
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
  struct mesh* m_out = new_mesh(elem_dim, mesh_get_rep(m), 0);
  refine_verts(m, m_out, src_dim, gen_offset_of_srcs);
  refine_ents(m, m_out, src_dim, gen_offset_of_srcs);
  loop_free(gen_offset_of_srcs);
  free_mesh(m);
  *p_m = m_out;
  return 1;
}
