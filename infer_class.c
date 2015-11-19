#include "infer_class.h"

#include "loop.h"
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
      if (dim == INVALID || vert_dim > dim) {
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
  *p_class_id_of_ents = id_of_ents;
}
