#include "swap_fit.hpp"

#include <cassert>
#include <cstdio>

#include "algebra.hpp"
#include "arrays.hpp"
#include "element_field.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "qr.hpp"
#include "size.hpp"

/* I tried to fit a linear polynomial here,
   but because these are tetrahedra around
   an edge, their centroids are nearly coplanar
   so the gradient along the edge is unstable.
   I'm switching to a dumb average, because
   diffusion is better than values outside
   the sensible range */

LOOP_KERNEL(swap_fit_cavity,
    unsigned ncomps,
    unsigned const* gen_offset_of_edges,
    unsigned const* elems_of_edges_offsets,
    unsigned const* elems_of_edges,
    double const* data_in,
    double* gen_data)

  if (gen_offset_of_edges[i] ==
      gen_offset_of_edges[i + 1])
    return;
  unsigned f = elems_of_edges_offsets[i];
  unsigned e = elems_of_edges_offsets[i + 1];
  for (unsigned j = 0; j < ncomps; ++j) {
    double avg = 0;
    for (unsigned k = f; k < e; ++k) {
      unsigned old_elem = elems_of_edges[k];
      avg += data_in[old_elem * ncomps + j];
    }
    avg /= (e - f);
    unsigned nf = gen_offset_of_edges[i];
    unsigned ne = gen_offset_of_edges[i + 1];
    for (unsigned gen_elem = nf; gen_elem < ne; ++gen_elem) {
      gen_data[gen_elem * ncomps + j] = avg;
    }
  }
}

static double* swap_fit_data(
    struct mesh* m,
    unsigned const* gen_offset_of_edges,
    struct const_tag* t)
{
  unsigned ncomps = t->ncomps;
  unsigned elem_dim = mesh_dim(m);
  unsigned nedges = mesh_count(m, 1);
  unsigned ngen_elems = array_at(gen_offset_of_edges, nedges);
  unsigned const* elems_of_edges_offsets =
    mesh_ask_up(m, 1, elem_dim)->offsets;
  unsigned const* elems_of_edges =
    mesh_ask_up(m, 1, elem_dim)->adj;
  double const* data_in = t->d.f64;
  double* gen_data = LOOP_MALLOC(double, ngen_elems * ncomps);
  LOOP_EXEC(swap_fit_cavity, nedges,
      ncomps,
      gen_offset_of_edges,
      elems_of_edges_offsets,
      elems_of_edges,
      data_in,
      gen_data);
  return gen_data;
}

static void swap_fit_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_edges,
    unsigned const* offset_of_same_elems,
    struct const_tag* t)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned nedges = mesh_count(m, 1);
  double* same_data = expand_array(t->d.f64, offset_of_same_elems,
      nelems, t->ncomps);
  unsigned nsame_elems = array_at(offset_of_same_elems, nelems);
  double* gen_data = swap_fit_data(m, gen_offset_of_edges, t);
  double* data_out = concat_arrays(same_data, gen_data,
      nsame_elems, array_at(gen_offset_of_edges, nedges), t->ncomps);
  loop_free(same_data);
  loop_free(gen_data);
  add_tag2(mesh_tags(m_out, elem_dim), TAG_F64, t->name, t->ncomps,
      t->transfer_type, data_out);
}

static int should_fit(struct mesh* m, unsigned tag_i)
{
  unsigned elem_dim = mesh_dim(m);
  struct const_tag* t = mesh_get_tag(m, elem_dim, tag_i);
  return t->type == TAG_F64 && t->transfer_type == OSH_TRANSFER_POINTWISE;
}

void swap_fit(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_edges,
    unsigned const* offset_of_same_elems)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned i;
  for (i = 0; i < mesh_count_tags(m, elem_dim); ++i)
    if (should_fit(m, i))
      break;
  if (i == mesh_count_tags(m, elem_dim))
    return;
  assert(elem_dim == 3);
  for (i = 0; i < mesh_count_tags(m, elem_dim); ++i) {
    struct const_tag* t = mesh_get_tag(m, elem_dim, i);
    if ((t->type == TAG_F64) && (t->transfer_type == OSH_TRANSFER_POINTWISE))
      swap_fit_tag(m, m_out, gen_offset_of_edges, offset_of_same_elems, t);
  }
}

