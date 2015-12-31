#include "refine_class.h"

#include <assert.h>

#include "arrays.h"
#include "infer_class.h"
#include "inherit.h"
#include "loop.h"
#include "mesh.h"

void refine_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned* offset_of_doms[4])
{
  if (!mesh_find_tag(m_in, 0, "class_dim"))
    return;
  if (mesh_get_rep(m_in) == MESH_REDUCED)
    mesh_ask_class(m_in, src_dim);
  for (unsigned prod_dim = 0; prod_dim <= mesh_dim(m_in); ++prod_dim) {
    if (mesh_get_rep(m_in) == MESH_REDUCED && prod_dim != 0)
      continue;
    unsigned ndoms[4];
    unsigned* prods_of_doms_offsets[4];
    setup_refine(m_in, src_dim, prod_dim, offset_of_doms,
        ndoms, prods_of_doms_offsets);
    inherit_uint_tag(m_in, m_out, prod_dim, ndoms, prods_of_doms_offsets,
        "class_dim");
    if (mesh_find_tag(m_in, 0, "class_id"))
      inherit_uint_tag(m_in, m_out, prod_dim, ndoms, prods_of_doms_offsets,
          "class_id");
    for (unsigned i = 0; i < 4; ++i)
      loop_free(prods_of_doms_offsets[i]);
  }
}
