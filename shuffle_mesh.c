#include "shuffle_mesh.h"

#include "arrays.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"

#if 0

static void shuffle_tags(struct mesh* in, struct mesh* out,
    unsigned dim, unsigned const* old_to_new)
{
  if (mesh_is_parallel(in))
    mesh_parallel_to_tags(in, dim);
  unsigned nents = mesh_count(in, dim);
  for (unsigned i = 0; i < mesh_count_tags(in, dim); ++i) {
    struct const_tag* t = mesh_get_tag(in, dim, i);
    void* vals_out = 0;
    switch (t->type) {
      case TAG_U8:
        break;
      case TAG_U32:
        vals_out = uints_shuffle(nents, t->d.u32, t->ncomps,
            old_to_new);
        break;
      case TAG_U64:
        vals_out = ulongs_shuffle(nents, t->d.u64, t->ncomps,
            old_to_new);
        break;
      case TAG_F64:
        vals_out = doubles_shuffle(nents, t->d.f64, t->ncomps,
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
  unsigned* verts_of_ents_out = uints_shuffle(nents, verts_of_ents,
      verts_per_ent, old_to_new_ents);
  LOOP_EXEC(remap_conn, nents * verts_per_ent, old_to_new_verts,
      verts_of_ents_out);
  shuffle_tags(m, m_out, ent_dim, old_to_new_ents);
}

#endif

LOOP_KERNEL(count_fan_ents,
    unsigned const* vert_num,
    unsigned const* ents_of_verts_offsets,
    unsigned const* ents_of_verts,
    unsigned const* verts_of_ents,
    unsigned verts_per_ent,
    unsigned* nfan_ents)
  unsigned n = vert_num[i];
  unsigned a = ents_of_verts_offsets[i];
  unsigned b = ents_of_verts_offsets[i + 1];
  unsigned nfans = 0;
  for (unsigned j = a; j < b; ++j) {
    unsigned ent = ents_of_verts[j];
    unsigned const* verts_of_ent = verts_of_ents + ent * verts_per_ent;
    unsigned k;
    for (k = 0; k < verts_per_ent; ++k) {
      unsigned oi = verts_of_ent[k];
      unsigned on = vert_num[oi];
      if (on < n)
        break;
    }
    if (k == verts_per_ent)
      ++nfans;
  }
  nfan_ents[i] = nfans;
}

LOOP_KERNEL(fill_fan_ents,
    unsigned const* vert_num,
    unsigned const* ents_of_verts_offsets,
    unsigned const* ents_of_verts,
    unsigned const* verts_of_ents,
    unsigned verts_per_ent,
    unsigned const* fan_ent_offsets,
    unsigned* fan_ents)
  unsigned n = vert_num[i];
  unsigned a = ents_of_verts_offsets[i];
  unsigned b = ents_of_verts_offsets[i + 1];
  unsigned* vert_fans = fan_ents + fan_ent_offsets[i];
  unsigned nfans = 0;
  for (unsigned j = a; j < b; ++j) {
    unsigned ent = ents_of_verts[j];
    unsigned const* verts_of_ent = verts_of_ents + ent * verts_per_ent;
    unsigned k;
    for (k = 0; k < verts_per_ent; ++k) {
      unsigned oi = verts_of_ent[k];
      unsigned on = vert_num[oi];
      if (on < n)
        break;
    }
    if (k == verts_per_ent)
      vert_fans[nfans++] = ent;
  }
}

LOOP_KERNEL(number_fan_ents,
    unsigned const* vert_num,
    unsigned const* fan_ent_offsets,
    unsigned const* fan_ents,
    unsigned const* new_vert_offsets,
    unsigned* ent_num)
  unsigned en = new_vert_offsets[vert_num[i]];
  unsigned a = fan_ent_offsets[i];
  unsigned b = fan_ent_offsets[i + 1];
  for (unsigned j = a; j < b; ++j) {
    unsigned ent = fan_ents[j];
    ent_num[ent] = en++;
  }
}

unsigned* number_ents(struct mesh* m,
    unsigned ent_dim, unsigned const* vert_num)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned nents = mesh_count(m, ent_dim);
  unsigned const* ents_of_verts_offsets = mesh_ask_up(m, 0, ent_dim)->offsets;
  unsigned const* ents_of_verts = mesh_ask_up(m, 0, ent_dim)->adj;
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned* nfan_ents = LOOP_MALLOC(unsigned, nverts);
  LOOP_EXEC(count_fan_ents, nverts, vert_num,
      ents_of_verts_offsets, ents_of_verts,
      verts_of_ents, verts_per_ent, nfan_ents);
  unsigned* fan_ent_offsets = uints_exscan(nfan_ents, nverts);
  assert(uints_at(fan_ent_offsets, nverts) == nents);
  unsigned* fan_ents = LOOP_MALLOC(unsigned, nents);
  LOOP_EXEC(fill_fan_ents, nverts, vert_num,
      ents_of_verts_offsets, ents_of_verts,
      verts_of_ents, verts_per_ent, fan_ent_offsets, fan_ents);
  unsigned* new_vert_nfans = uints_shuffle(nverts, nfan_ents, 1, vert_num);
  loop_free(nfan_ents);
  unsigned* new_vert_offsets = uints_exscan(new_vert_nfans, nverts);
  loop_free(new_vert_nfans);
  assert(uints_at(new_vert_offsets, nverts) == nents);
  unsigned* ent_num = LOOP_MALLOC(unsigned, nents);
  LOOP_EXEC(number_fan_ents, nverts, vert_num, fan_ent_offsets, fan_ents,
      new_vert_offsets, ent_num);
  loop_free(fan_ents);
  loop_free(fan_ent_offsets);
  loop_free(new_vert_offsets);
  return ent_num;
}
