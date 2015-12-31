#include "refine_class.h"

#include "arrays.h"
#include "infer_class.h"
#include "inherit.h"
#include "mesh.h"

void refine_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4])
{
  if (!mesh_find_tag(m_in, 0, "class_dim"))
    return;
  if (mesh_get_rep(m_in) == MESH_REDUCED) {
    if (prod_dim != 0)
      return;
    mesh_ask_class(m_in, src_dim);
  }
  inherit_uint_tag(m_in, m_out, prod_dim, ndoms, prods_of_doms_offsets,
      "class_dim");
  if (mesh_find_tag(m_in, 0, "class_id"))
    inherit_uint_tag(m_in, m_out, prod_dim, ndoms, prods_of_doms_offsets,
        "class_id");
}
