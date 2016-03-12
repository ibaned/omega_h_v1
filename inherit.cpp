#include "inherit.hpp"

#include <cassert>

#include "arrays.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "refine_topology.hpp"
#include "tables.hpp"

void make_ngen_from_doms(
    unsigned const ndoms[4],
    unsigned* prods_of_doms_offsets[4],
    unsigned ngen[4])
{
  for (unsigned i = 0; i < 4; ++i) {
    if (ndoms[i] && prods_of_doms_offsets[i])
      ngen[i] = uints_at(prods_of_doms_offsets[i], ndoms[i]);
    else
      ngen[i] = 0;
  }
}

void make_ngen_offsets(
    unsigned const ngen[4],
    /* out: */
    unsigned ngen_offsets[5])
{
  ngen_offsets[0] = 0;
  for (unsigned i = 0; i < 4; ++i)
    ngen_offsets[i + 1] = ngen_offsets[i] + ngen[i];
}

#define GENERIC_CONCAT_INHERITED(T, name) \
T* concat_##name##_inherited( \
    unsigned width, \
    unsigned const ngen_offsets[5], \
    T* gen_data[4]) \
{ \
  T* out_data = LOOP_MALLOC(T, ngen_offsets[4] * width); \
  for (unsigned i = 0; i < 4; ++i) { \
    array_memcpy<T>( \
        out_data + ngen_offsets[i] * width, \
        gen_data[i], \
        (ngen_offsets[i + 1] - ngen_offsets[i]) * width); \
    loop_free(gen_data[i]); \
  } \
  return out_data; \
}

GENERIC_CONCAT_INHERITED(unsigned, uints)
GENERIC_CONCAT_INHERITED(double, doubles)

static void inherit_uints_direct(
    unsigned width,
    unsigned ndoms[4],
    unsigned const* dom_data[4],
    unsigned* prods_of_doms_offsets[4],
    /* out: */
    unsigned* gen_data[4])
{
  for (unsigned i = 0; i < 4; ++i) {
    if (ndoms[i] && dom_data[i] && prods_of_doms_offsets[i]) {
      gen_data[i] = expand_array(ndoms[i], width,
          dom_data[i], prods_of_doms_offsets[i]);
    } else {
      gen_data[i] = 0;
    }
  }
}

static void set_zero(
    unsigned dom_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4])
{
  ndoms[dom_dim] = 0;
  prods_of_doms_offsets[dom_dim] = 0;
}

static void set_expanded_offsets(
    struct mesh* m,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned dom_dim,
    unsigned* gen_dom_offsets[4],
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4])
{
  ndoms[dom_dim] = mesh_count(m, dom_dim);
  unsigned prods_per_dom = get_prods_per_dom(dom_dim, src_dim, prod_dim);
  prods_of_doms_offsets[dom_dim] = uints_scale(gen_dom_offsets[dom_dim],
      ndoms[dom_dim] + 1, prods_per_dom);
}

void setup_refine(
    struct mesh* m,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned* gen_dom_offsets[4],
    /* out: */
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4],
    unsigned ngen_offsets[5])
{
  unsigned nprods = mesh_count(m, prod_dim);
  ndoms[0] = nprods;
  prods_of_doms_offsets[0] = uints_negate_offsets(
      gen_dom_offsets[prod_dim], nprods);
  if (mesh_get_rep(m) == MESH_FULL) {
    unsigned start = (src_dim > prod_dim) ? src_dim : prod_dim;
    for (unsigned d = 1; d < start; ++d)
      set_zero(d, ndoms, prods_of_doms_offsets);
    for (unsigned d = start; d <= mesh_dim(m); ++d)
      set_expanded_offsets(m, src_dim, prod_dim, d, gen_dom_offsets,
          ndoms, prods_of_doms_offsets);
  } else {
    unsigned dom_dim = (prod_dim == 0) ? src_dim : mesh_dim(m);
    for (unsigned d = 1; d <= mesh_dim(m); ++d)
      if (d != dom_dim)
        set_zero(d, ndoms, prods_of_doms_offsets);
    set_expanded_offsets(m, src_dim, prod_dim, dom_dim, gen_dom_offsets,
        ndoms, prods_of_doms_offsets);
  }
  for (unsigned d = mesh_dim(m) + 1; d < 4; ++d)
    set_zero(d, ndoms, prods_of_doms_offsets);
  unsigned ngen[4];
  make_ngen_from_doms(ndoms, prods_of_doms_offsets, ngen);
  make_ngen_offsets(ngen, ngen_offsets);
}

void setup_coarsen(
    struct mesh* m,
    unsigned ent_dim,
    unsigned* gen_offset_of_ents,
    unsigned* offset_of_same_ents,
    /* out: */
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4])
{
  unsigned nents = mesh_count(m, ent_dim);
  ndoms[0] = nents;
  prods_of_doms_offsets[0] = offset_of_same_ents;
  for (unsigned d = 1; d < 4; ++d)
    if (d != ent_dim)
      set_zero(d, ndoms, prods_of_doms_offsets);
  ndoms[ent_dim] = nents;
  prods_of_doms_offsets[ent_dim] = gen_offset_of_ents;
}

void setup_swap(
    struct mesh* m,
    unsigned ent_dim,
    unsigned* gen_offset_of_edges,
    unsigned* offset_of_same_ents,
    /* out: */
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4])
{
  ndoms[0] = mesh_count(m, ent_dim);
  prods_of_doms_offsets[0] = offset_of_same_ents;
  ndoms[1] = mesh_count(m, 1);
  prods_of_doms_offsets[1] = gen_offset_of_edges;
  for (unsigned d = 2; d < 4; ++d)
    set_zero(d, ndoms, prods_of_doms_offsets);
}

static struct const_tag* setup_uint_tag(
    struct mesh* m,
    char const* name,
    unsigned prod_dim,
    unsigned const* dom_data[4])
{
  struct const_tag* t = mesh_find_tag(m, prod_dim, name);
  dom_data[0] = t->d.u32;
  for (unsigned d = 1; d <= mesh_dim(m); ++d) {
    if (mesh_find_tag(m, d, name))
      dom_data[d] = mesh_find_tag(m, d, name)->d.u32;
    else
      dom_data[d] = 0;
  }
  for (unsigned d = mesh_dim(m) + 1; d < 4; ++d)
    dom_data[d] = 0;
  return t;
}

static void inherit_uint_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned prod_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4],
    char const* name)
{
  unsigned const* dom_data[4];
  struct const_tag* t;
  t = setup_uint_tag(m, name, prod_dim, dom_data);
  unsigned ngen[4];
  make_ngen_from_doms(ndoms, prods_of_doms_offsets, ngen);
  unsigned ngen_offsets[5];
  make_ngen_offsets(ngen, ngen_offsets);
  unsigned* gen_data[4];
  inherit_uints_direct(t->ncomps, ndoms, dom_data, prods_of_doms_offsets,
      gen_data);
  unsigned* data_out = concat_uints_inherited(t->ncomps, ngen_offsets, gen_data);
  mesh_add_tag(m_out, prod_dim, TAG_U32, t->name, t->ncomps, data_out);
}

void inherit_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned prod_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4])
{
  if (!mesh_find_tag(m_in, 0, "class_dim"))
    return;
  assert(mesh_get_rep(m_in) == MESH_FULL);
  inherit_uint_tag(m_in, m_out, prod_dim, ndoms, prods_of_doms_offsets,
      "class_dim");
  if (mesh_find_tag(m_in, 0, "class_id"))
    inherit_uint_tag(m_in, m_out, prod_dim, ndoms, prods_of_doms_offsets,
        "class_id");
}

void concat_verts_of_ents(
    unsigned ent_dim,
    unsigned nents,
    unsigned ngen_ents,
    unsigned const* verts_of_ents,
    unsigned const* offset_of_same_ents,
    unsigned const* verts_of_gen_ents,
    unsigned* p_nents_out,
    unsigned** p_verts_of_ents_out)
{
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned nsame_ents = uints_at(offset_of_same_ents, nents);
  unsigned nents_out = nsame_ents + ngen_ents;
  unsigned* verts_of_same_ents = expand_array(nents, verts_per_ent,
      verts_of_ents, offset_of_same_ents);
  unsigned* verts_of_ents_out = concat_arrays(verts_per_ent,
      verts_of_same_ents, nsame_ents,
      verts_of_gen_ents, ngen_ents);
  loop_free(verts_of_same_ents);
  *p_nents_out = nents_out;
  *p_verts_of_ents_out = verts_of_ents_out;
}
