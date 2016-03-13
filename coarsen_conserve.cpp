#include "coarsen_conserve.hpp"

#include <cassert>

#include "algebra.hpp"
#include "arrays.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "size.hpp"

#define MAX_NCOMPS 32

LOOP_KERNEL(coarsen_conserve_cavity,
    unsigned ncomps,
    unsigned nsame_elems,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_offset_of_elems,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    double const* data_in,
    double const* new_elem_sizes,
    double* gen_data)

  if (gen_offset_of_verts[i] ==
      gen_offset_of_verts[i + 1])
    return;
  unsigned f = elems_of_verts_offsets[i];
  unsigned e = elems_of_verts_offsets[i + 1];
  double sum[MAX_NCOMPS] = {0};
  for (unsigned j = f; j < e; ++j) {
    unsigned elem = elems_of_verts[j];
    add_vectors(sum, data_in + elem * ncomps, sum, ncomps);
  }
  double new_cavity_size = 0;
  for (unsigned j = f; j < e; ++j) {
    unsigned elem = elems_of_verts[j];
    if (gen_offset_of_elems[elem] ==
        gen_offset_of_elems[elem + 1])
      continue;
    unsigned new_elem = gen_offset_of_elems[elem] + nsame_elems;
    new_cavity_size += new_elem_sizes[new_elem];
  }
  for (unsigned j = f; j < e; ++j) {
    unsigned elem = elems_of_verts[j];
    if (gen_offset_of_elems[elem] ==
        gen_offset_of_elems[elem + 1])
      continue;
    unsigned gen_elem = gen_offset_of_elems[elem];
    unsigned new_elem = gen_elem + nsame_elems;
    scale_vector(sum, new_elem_sizes[new_elem] / new_cavity_size,
        gen_data + gen_elem * ncomps, ncomps);
  }
}

static double* coarsen_conserve_data(
    struct mesh* m,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_offset_of_elems,
    unsigned nsame_elems,
    double const* new_elem_sizes,
    struct const_tag* t)
{
  unsigned ncomps = t->ncomps;
  assert(ncomps <= MAX_NCOMPS);
  unsigned elem_dim = mesh_dim(m);
  unsigned nverts = mesh_count(m, 0);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* elems_of_verts_offsets =
    mesh_ask_up(m, 0, elem_dim)->offsets;
  unsigned const* elems_of_verts =
    mesh_ask_up(m, 0, elem_dim)->adj;
  double const* data_in = t->d.f64;
  unsigned ngen_elems = array_at(gen_offset_of_elems, nelems);
  double* gen_data = LOOP_MALLOC(double, ngen_elems * ncomps);
  LOOP_EXEC(coarsen_conserve_cavity, nverts,
      ncomps,
      nsame_elems,
      gen_offset_of_verts,
      gen_offset_of_elems,
      elems_of_verts_offsets,
      elems_of_verts,
      data_in,
      new_elem_sizes,
      gen_data);
  return gen_data;
}

static void coarsen_conserve_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_offset_of_elems,
    unsigned const* offset_of_same_elems,
    double const* new_elem_sizes,
    struct const_tag* t)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  double* same_data = expand_array(t->d.f64, offset_of_same_elems,
      nelems, t->ncomps);
  unsigned nsame_elems = array_at(offset_of_same_elems, nelems);
  double* gen_data = coarsen_conserve_data(m, gen_offset_of_verts,
      gen_offset_of_elems, nsame_elems, new_elem_sizes, t);
  double* data_out = concat_arrays(t->ncomps,
      same_data, nsame_elems,
      gen_data, array_at(gen_offset_of_elems, nelems));
  loop_free(same_data);
  loop_free(gen_data);
  add_tag2(mesh_tags(m_out, elem_dim), TAG_F64, t->name, t->ncomps,
      t->transfer_type, data_out);
}

void coarsen_conserve(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_offset_of_elems,
    unsigned const* offset_of_same_elems)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned i;
  for (i = 0; i < mesh_count_tags(m, elem_dim); ++i)
    if (mesh_get_tag(m, elem_dim, i)->type == TAG_F64)
      break;
  if (i == mesh_count_tags(m, elem_dim))
    return;
  double* new_elem_sizes = mesh_element_sizes(m_out);
  for (i = 0; i < mesh_count_tags(m, elem_dim); ++i) {
    struct const_tag* t = mesh_get_tag(m, elem_dim, i);
    if ((t->type == TAG_F64) && (t->transfer_type == OSH_TRANSFER_CONSERVE))
      coarsen_conserve_tag(m, m_out, gen_offset_of_verts,
          gen_offset_of_elems, offset_of_same_elems,
          new_elem_sizes, t);
  }
  loop_free(new_elem_sizes);
}
