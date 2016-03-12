#include "shuffle_mesh.hpp"

#include "arrays.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "reorder.hpp"
#include "tables.hpp"

static void shuffle_tags(struct mesh* in, struct mesh* out,
    unsigned dim, unsigned const* old_to_new)
{
  if (mesh_is_parallel(in))
    mesh_tag_globals(in, dim);
  unsigned nents = mesh_count(in, dim);
  for (unsigned i = 0; i < mesh_count_tags(in, dim); ++i) {
    struct const_tag* t = mesh_get_tag(in, dim, i);
    void* vals_out = 0;
    switch (t->type) {
      case TAG_U8:
        break;
      case TAG_U32:
        vals_out = shuffle_array(nents, t->d.u32, t->ncomps,
            old_to_new);
        break;
      case TAG_U64:
        vals_out = shuffle_array(nents, t->d.u64, t->ncomps,
            old_to_new);
        break;
      case TAG_F64:
        vals_out = shuffle_array(nents, t->d.f64, t->ncomps,
            old_to_new);
        break;
    }
    mesh_add_tag(out, dim, t->type, t->name, t->ncomps, vals_out);
  }
  if (mesh_is_parallel(in)) {
    mesh_parallel_untag(in, dim);
    mesh_parallel_from_tags(out, dim);
  }
}

/* FIXME: duplicated in coarsen_common.c */
LOOP_KERNEL(remap_conn,
    unsigned const* old_to_new_verts,
    unsigned* verts_of_ents_out)
  verts_of_ents_out[i] = old_to_new_verts[verts_of_ents_out[i]];
}

static void shuffle_ents(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ent_dim,
    unsigned const* old_to_new_ents,
    unsigned const* old_to_new_verts)
{
  unsigned nents = mesh_count(m, ent_dim);
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned* verts_of_ents_out = shuffle_array(nents, verts_of_ents,
      verts_per_ent, old_to_new_ents);
  LOOP_EXEC(remap_conn, nents * verts_per_ent, old_to_new_verts,
      verts_of_ents_out);
  mesh_set_ents(m_out, ent_dim, nents, verts_of_ents_out);
  shuffle_tags(m, m_out, ent_dim, old_to_new_ents);
}

void shuffle_mesh(struct mesh* m, unsigned const* old_to_new_verts)
{
  unsigned dim = mesh_dim(m);
  struct mesh* m_out = new_mesh(dim, mesh_get_rep(m), mesh_is_parallel(m));
  mesh_set_ents(m_out, 0, mesh_count(m, 0), 0);
  shuffle_tags(m, m_out, 0, old_to_new_verts);
  for (unsigned d = 1; d <= dim; ++d) {
    if (!mesh_has_dim(m, d))
      continue;
    unsigned* old_to_new_ents = number_ents(m, d, old_to_new_verts);
    shuffle_ents(m, m_out, d, old_to_new_ents, old_to_new_verts);
    loop_free(old_to_new_ents);
  }
  overwrite_mesh(m, m_out);
}
