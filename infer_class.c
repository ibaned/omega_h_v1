#include "infer_class.h"

#include <assert.h>

#include "loop.h"
#include "mesh.h"
#include "tables.h"

void infer_class(
    unsigned ent_dim,
    unsigned nents,
    unsigned const* verts_of_ents,
    unsigned const* class_dim_of_verts,
    unsigned const* class_id_of_verts,
    unsigned** p_class_dim_of_ents,
    unsigned** p_class_id_of_ents)
{
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned* dim_of_ents = LOOP_MALLOC(unsigned, nents);
  unsigned* id_of_ents = 0;
  if (class_id_of_verts)
    id_of_ents = LOOP_MALLOC(unsigned, nents);
  for (unsigned i = 0; i < nents; ++i) {
    unsigned const* verts_of_ent = verts_of_ents + i * verts_per_ent;
    unsigned dim = INVALID;
    unsigned id = INVALID;
    for (unsigned j = 0; j < verts_per_ent; ++j) {
      unsigned vert = verts_of_ent[j];
      unsigned vert_dim = class_dim_of_verts[vert];
      if (j == 0 || vert_dim > dim) {
        dim = vert_dim;
        if (class_id_of_verts)
          id = class_id_of_verts[vert];
      }
    }
    dim_of_ents[i] = dim;
    if (class_id_of_verts)
      id_of_ents[i] = id;
  }
  *p_class_dim_of_ents = dim_of_ents;
  if (class_id_of_verts)
    *p_class_id_of_ents = id_of_ents;
}

void mesh_ask_class(struct mesh* m, unsigned dim)
{
  if (mesh_find_tag(m, dim, "class_dim"))
    return;
  assert(mesh_find_tag(m, 0, "class_dim"));
  unsigned nents = mesh_count(m, dim);
  unsigned const* verts_of_ents = mesh_ask_down(m, dim, 0);
  unsigned const* class_dim_of_verts = mesh_find_tag(m, 0, "class_dim")->d.u32;
  unsigned const* class_id_of_verts = 0;
  if (mesh_find_tag(m, 0, "class_id"))
    class_id_of_verts = mesh_find_tag(m, 0, "class_id")->d.u32;
  unsigned* class_dim_of_ents;
  unsigned* class_id_of_ents;
  infer_class(dim, nents, verts_of_ents, class_dim_of_verts, class_id_of_verts,
      &class_dim_of_ents, &class_id_of_ents);
  mesh_add_tag(m, dim, TAG_U32, "class_dim", 1, class_dim_of_ents);
  if (class_id_of_ents)
    mesh_add_tag(m, dim, TAG_U32, "class_id", 1, class_id_of_ents);
}

unsigned const* mesh_ask_class_dim(struct mesh* m, unsigned dim)
{
  mesh_ask_class(m, dim);
  return mesh_find_tag(m, dim, "class_dim")->d.u32;
}

unsigned const* mesh_ask_class_id(struct mesh* m, unsigned dim)
{
  mesh_ask_class(m, dim);
  if (!mesh_find_tag(m, dim, "class_id"))
    return 0;
  return mesh_find_tag(m, dim, "class_id")->d.u32;
}
