#include "refine_class.h"

#include <assert.h>

#include "arrays.h"
#include "infer_class.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "subset.h"

static void inherit_uint_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned* offset_of_doms[4],
    unsigned ngen_ents[4][4],
    unsigned gen_offsets[4][5],
    char const* name)
{
  unsigned elem_dim = mesh_dim(m);
  for (unsigned prod_dim = 0; prod_dim <= elem_dim; ++prod_dim) {
    if ((mesh_get_rep(m) == MESH_REDUCED) &&
        !((prod_dim == 0)||(prod_dim == elem_dim)))
      continue;
    unsigned nprods = gen_offsets[prod_dim][4];
    struct const_tag* prod_tag = mesh_find_tag(m, prod_dim, name);
    if (!prod_tag)
      continue;
    unsigned* prod_data = LOOP_MALLOC(unsigned, nprods * prod_tag->ncomps);
    unsigned nold = mesh_count(m, prod_dim);
    unsigned* offset_of_same = uints_negate_offsets(offset_of_doms[prod_dim],
        nold);
    assert(offset_of_same[nold] == gen_offsets[prod_dim][0]);
    uints_subset_into(nold, prod_tag->ncomps,
        prod_tag->d.u32, offset_of_same, prod_data);
    loop_free(offset_of_same);
    for (unsigned dom_dim = prod_dim; dom_dim <= elem_dim; ++dom_dim) {
      if (dom_dim < src_dim)
        continue;
      if (mesh_get_rep(m) == MESH_REDUCED && dom_dim != elem_dim)
        continue;
      unsigned ndoms = mesh_count(m, dom_dim);
      unsigned nsplit_doms = offset_of_doms[dom_dim][ndoms];
      assert(ngen_ents[prod_dim][dom_dim] % nsplit_doms == 0);
      unsigned prods_per_dom = ngen_ents[prod_dim][dom_dim] / nsplit_doms;
      if (!prods_per_dom)
        continue;
      struct const_tag* dom_tag = mesh_find_tag(m, dom_dim, name);
      assert(dom_tag->ncomps == prod_tag->ncomps);
      assert(dom_tag->type == prod_tag->type);
      unsigned* offsets = LOOP_MALLOC(unsigned, ndoms + 1);
      for (unsigned i = 0; i <= ndoms; ++i)
        offsets[i] = (offset_of_doms[dom_dim][i] * prods_per_dom) +
          gen_offsets[prod_dim][dom_dim];
      uints_expand_into(ndoms, dom_tag->d.u32, dom_tag->ncomps,
          offsets, prod_data);
      loop_free(offsets);
    }
    mesh_add_tag(m_out, prod_dim, prod_tag->type, prod_tag->name,
        prod_tag->ncomps, prod_data);
  }
}

void refine_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned* offset_of_doms[4],
    unsigned ngen_ents[4][4],
    unsigned gen_offsets[4][5])
{
  if (mesh_get_rep(m_in) == MESH_REDUCED &&
      mesh_find_tag(m_in, 0, "class_dim"))
    mesh_ask_class(m_in, src_dim);
  inherit_uint_tag(m_in, m_out, src_dim, offset_of_doms, ngen_ents,
      gen_offsets, "class_dim");
  inherit_uint_tag(m_in, m_out, src_dim, offset_of_doms, ngen_ents,
      gen_offsets, "class_id");
}
