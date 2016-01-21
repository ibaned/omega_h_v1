#include "subset.h"

#include <string.h>

#include "arrays.h"
#include "ints.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"

void tags_subset(struct mesh* in, struct mesh* out,
    unsigned dim, unsigned const* offsets)
{
  unsigned nverts = mesh_count(in, dim);
  for (unsigned i = 0; i < mesh_count_tags(in, dim); ++i) {
    struct const_tag* t = mesh_get_tag(in, dim, i);
    void* vals_out = 0;
    switch (t->type) {
      case TAG_U8:
        vals_out = uchars_expand(nverts, t->ncomps, t->d.u8,
            offsets);
        break;
      case TAG_U32:
        vals_out = uints_expand(nverts, t->ncomps, t->d.u32,
            offsets);
        break;
      case TAG_U64:
        vals_out = ulongs_expand(nverts, t->ncomps, t->d.u64,
            offsets);
        break;
      case TAG_F64:
        vals_out = doubles_expand(nverts, t->ncomps, t->d.f64,
            offsets);
        break;
    }
    mesh_add_tag(out, dim, t->type, t->name, t->ncomps, vals_out);
  }
  if (mesh_is_parallel(in)) {
    mesh_parallel_untag(in, dim);
    mesh_parallel_from_tags(out, dim);
  }
}

static void subset_ents(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ent_dim,
    unsigned const* ent_offsets,
    unsigned const* vert_offsets)
{
  unsigned nents = mesh_count(m, ent_dim);
  unsigned nents_out = ent_offsets[nents];
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned* verts_of_ents_out = uints_expand(nents, verts_per_ent,
      verts_of_ents, ent_offsets);
  for (unsigned i = 0; i < nents_out * verts_per_ent; ++i)
    verts_of_ents_out[i] = vert_offsets[verts_of_ents_out[i]];
  mesh_set_ents(m_out, ent_dim, nents_out, verts_of_ents_out);
  tags_subset(m, m_out, ent_dim, ent_offsets);
}

struct mesh* subset_mesh(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* offsets)
{
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned* marked_elems = uints_unscan(offsets, nelems);
  unsigned* to_free[4] = {0};
  unsigned const* ent_offsets[4] = {0};
  ent_offsets[elem_dim] = offsets;
  for (unsigned d = 0; d < elem_dim; ++d) {
    if (!mesh_has_dim(m, d))
      continue;
    unsigned nents = mesh_count(m, d);
    unsigned* marked_ents = mesh_mark_down(m, elem_dim, d, marked_elems);
    ent_offsets[d] = to_free[d] = uints_exscan(marked_ents, nents);
    loop_free(marked_ents);
  }
  loop_free(marked_elems);
  unsigned nverts = mesh_count(m, 0);
  unsigned nverts_out = ent_offsets[0][nverts];
  struct mesh* m_out = new_mesh(elem_dim, mesh_get_rep(m), mesh_is_parallel(m));
  mesh_set_ents(m_out, 0, nverts_out, 0);
  tags_subset(m, m_out, 0, ent_offsets[0]);
  for (unsigned d = 1; d <= elem_dim; ++d)
    if (mesh_has_dim(m, d))
      subset_ents(m, m_out, d, ent_offsets[d], ent_offsets[0]);
  for (unsigned d = 0; d < elem_dim; ++d)
    loop_free(to_free[d]);
  return m_out;
}

void subset_verts_of_doms(
    struct mesh* m,
    unsigned dom_dim,
    unsigned const* offset_of_doms,
    unsigned* verts_of_prods)
{
  unsigned const* verts_of_doms = mesh_ask_down(m, dom_dim, 0);
  unsigned* offset_of_same = uints_negate_offsets(offset_of_doms,
      mesh_count(m, dom_dim));
  uints_expand_into(mesh_count(m, dom_dim), the_down_degrees[dom_dim][0],
      verts_of_doms, offset_of_same, verts_of_prods);
  loop_free(offset_of_same);
}
