#include "refine_fit.hpp"

#include "arrays.hpp"
#include "inherit.hpp"
#include "loop.hpp"
#include "mesh.hpp"

LOOP_KERNEL(refine_fit_elem,
    unsigned const* prods_of_doms_offsets,
    unsigned width,
    double const* data_in,
    double* data_out)
  unsigned f = prods_of_doms_offsets[i];
  unsigned e = prods_of_doms_offsets[i + 1];
  for (unsigned j = f; j < e; ++j)
    for (unsigned k = 0; k < width; ++k)
      data_out[j * width + k] = data_in[i * width + k];
}


static double* refine_fit_data(
    unsigned nelems,
    unsigned const* prods_of_doms_offsets,
    struct const_tag* t)
{
  unsigned width = t->ncomps;
  double const* data_in = t->d.f64;
  unsigned ngen_elems = array_at(prods_of_doms_offsets, nelems);
  double* data_out = LOOP_MALLOC(double, width * ngen_elems);
  LOOP_EXEC(refine_fit_elem, nelems,
      prods_of_doms_offsets, width, data_in, data_out);
  return data_out;
}

static void refine_fit_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4],
    struct const_tag* t)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned ngen[4];
  make_ngen_from_doms(ndoms, prods_of_doms_offsets, ngen);
  unsigned ngen_offsets[5];
  make_ngen_offsets(ngen, ngen_offsets);
  double* gen_data[4] = {0};
  gen_data[0] = expand_array(ndoms[0], t->ncomps, t->d.f64,
      prods_of_doms_offsets[0]);
  gen_data[elem_dim] = refine_fit_data(mesh_count(m, elem_dim),
      prods_of_doms_offsets[elem_dim], t);
  double* data_out = concat_inherited(t->ncomps, ngen_offsets, gen_data);
  add_tag2(mesh_tags(m_out, elem_dim), TAG_F64, t->name, t->ncomps,
      t->transfer_type, data_out);
}

void refine_fit(
    struct mesh* m,
    struct mesh* m_out,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4])
{
  unsigned elem_dim = mesh_dim(m);
  for (unsigned i = 0; i < mesh_count_tags(m, elem_dim); ++i) {
    struct const_tag* t = mesh_get_tag(m, elem_dim, i);
    if ((t->type == TAG_F64) && (t->transfer_type == OSH_TRANSFER_POINTWISE))
      refine_fit_tag(m, m_out, ndoms, prods_of_doms_offsets, t);
  }
}
