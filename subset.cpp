#include "subset.hpp"

#include <cstring>

#include "arrays.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mark.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "tables.hpp"

namespace omega_h {

void tags_subset(struct mesh* in, struct mesh* out,
    unsigned dim, unsigned const* offsets)
{
  unsigned nents = mesh_count(in, dim);
  for (unsigned i = 0; i < mesh_count_tags(in, dim); ++i) {
    struct const_tag* t = mesh_get_tag(in, dim, i);
    void* vals_out = 0;
    switch (t->type) {
      case TAG_U8:
        vals_out = expand_array(t->d.u8, offsets, nents, t->ncomps);
        break;
      case TAG_U32:
        vals_out = expand_array(t->d.u32, offsets, nents, t->ncomps);
        break;
      case TAG_U64:
        vals_out = expand_array(t->d.u64, offsets, nents, t->ncomps);
        break;
      case TAG_F64:
        vals_out = expand_array(t->d.f64, offsets, nents, t->ncomps);
        break;
    }
    add_tag2(mesh_tags(out, dim), t->type, t->name, t->ncomps,
        t->transfer_type, vals_out);
  }
  if (mesh_is_parallel(in)) {
    mesh_parallel_untag(in, dim);
    mesh_parallel_from_tags(out, dim);
  }
}

LOOP_KERNEL(remap_conn,
    unsigned const* offset_of_same_verts,
    unsigned* verts_of_ents_out)
  verts_of_ents_out[i] = offset_of_same_verts[verts_of_ents_out[i]];
}

static void subset_ents(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ent_dim,
    unsigned const* ent_offsets,
    unsigned const* vert_offsets)
{
  unsigned nents = mesh_count(m, ent_dim);
  unsigned nents_out = array_at(ent_offsets, nents);
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned* verts_of_ents_out = expand_array(verts_of_ents,
      ent_offsets, nents, verts_per_ent);
  LOOP_EXEC(remap_conn, nents_out * verts_per_ent, vert_offsets,
      verts_of_ents_out);
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
    unsigned* marked_ents = mesh_mark_down_local(m, elem_dim, d, marked_elems);
    ent_offsets[d] = to_free[d] = uints_exscan(marked_ents, nents);
    loop_free(marked_ents);
  }
  loop_free(marked_elems);
  unsigned nverts = mesh_count(m, 0);
  unsigned nverts_out = array_at(ent_offsets[0], nverts);
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
  expand_into(verts_of_prods, verts_of_doms, offset_of_same,
      mesh_count(m, dom_dim), the_down_degrees[dom_dim][0]);
  loop_free(offset_of_same);
}

}
