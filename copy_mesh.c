#include "copy_mesh.h"

#include "arrays.h"
#include "ints.h"
#include "doubles.h"
#include "mesh.h"
#include "tables.h"

void copy_tags(struct tags* a, struct tags* b, unsigned n)
{
  for (unsigned i = 0; i < count_tags(a); ++i) {
    struct const_tag* t = get_tag(a, i);
    void* data = 0;
    switch (t->type) {
      case TAG_U8:  data = uchars_copy(t->d.u8, n * t->ncomps);
                    break;
      case TAG_U32: data = uints_copy(t->d.u32, n * t->ncomps);
                    break;
      case TAG_U64: data = ulongs_copy(t->d.u64, n * t->ncomps);
                    break;
      case TAG_F64: data = doubles_copy(t->d.f64, n * t->ncomps);
                    break;
    };
    add_tag(b, t->type, t->name, t->ncomps, data);
  }
}

struct mesh* copy_mesh(struct mesh* a)
{
  unsigned dim = mesh_dim(a);
  struct mesh* b = new_mesh(dim, mesh_get_rep(a));
  for (unsigned d = 0; d <= dim; ++d) {
    if (!mesh_has_dim(a, d))
      continue;
    unsigned nents = mesh_count(a, d);
    unsigned verts_per_ent = the_down_degrees[d][0];
    unsigned* verts_of_ents = 0;
    if (d)
      verts_of_ents = uints_copy(mesh_ask_down(a, d, 0), nents * verts_per_ent);
    mesh_set_ents(b, d, nents, verts_of_ents);
    struct tags* tsa = mesh_tags(a, d);
    struct tags* tsb = mesh_tags(b, d);
    copy_tags(tsa, tsb, nents);
  }
  return b;
}
