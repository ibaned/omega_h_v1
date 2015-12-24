#include "refine_class.h"

#include <assert.h>

#include "arrays.h"
#include "infer_class.h"
#include "loop.h"
#include "mesh.h"
#include "subset.h"

void refine_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned const* gen_offset_of_srcs)
{
  /* if the input vertices are not classified there is no hope of
     classifying the output vertices, since most of them are the same */
  if (!mesh_find_tag(m_in, 0, "class_dim"))
    return;
  assert(mesh_find_tag(m_in, 0, "class_id"));
  unsigned const* class_dim_of_verts = mesh_find_tag(
      m_in, 0, "class_dim")->d.u32;
  unsigned const* class_id_of_verts = mesh_find_tag(
      m_in, 0, "class_id")->d.u32;
  /* either use infer_class to compute the classifications of the
     source-dimension entities or accept a pre-existing classification
     tag on said entities */
  unsigned const* class_dim_of_srcs = mesh_ask_class_dim(m_in, src_dim);
  unsigned const* class_id_of_srcs = mesh_ask_class_id(m_in, src_dim);
  unsigned nsrcs = mesh_count(m_in, src_dim);
  /* get the subset of the source-dimension entities that generated
     new vertices */
  unsigned* gen_class_dim = uints_subset(nsrcs, 1, class_dim_of_srcs,
      gen_offset_of_srcs);
  unsigned* gen_class_id = uints_subset(nsrcs, 1, class_id_of_srcs,
      gen_offset_of_srcs);
  unsigned nverts = mesh_count(m_in, 0);
  unsigned ngen = gen_offset_of_srcs[nsrcs];
  /* concatenate the new vertices with the existing ones */
  unsigned* class_dim_out = concat_uints(1, class_dim_of_verts, nverts,
      gen_class_dim, ngen);
  unsigned* class_id_out = concat_uints(1, class_id_of_verts, nverts,
      gen_class_id, ngen);
  loop_free(gen_class_dim);
  loop_free(gen_class_id);
  /* put the result on the output mesh */
  mesh_add_tag(m_out, 0, TAG_U32, "class_dim", 1, class_dim_out);
  mesh_add_tag(m_out, 0, TAG_U32, "class_id", 1, class_id_out);
}
