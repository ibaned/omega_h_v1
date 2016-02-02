#include "swap_conserve.h"

#include <assert.h>

#include "algebra.h"
#include "arrays.h"
#include "loop.h"
#include "mesh.h"
#include "size.h"

#define MAX_NCOMPS 32

static double* swap_conserve_data(
    struct mesh* m,
    unsigned const* gen_offset_of_edges,
    unsigned nsame_elems,
    double const* new_elem_sizes,
    struct const_tag* t)
{
  unsigned ncomps = t->ncomps;
  assert(ncomps <= MAX_NCOMPS);
  unsigned elem_dim = mesh_dim(m);
  unsigned nedges = mesh_count(m, 1);
  unsigned ngen_elems = uints_at(gen_offset_of_edges, nedges);
  unsigned const* elems_of_edges_offsets =
    mesh_ask_up(m, 1, elem_dim)->offsets;
  unsigned const* elems_of_edges =
    mesh_ask_up(m, 1, elem_dim)->adj;
  double const* data_in = t->d.f64;
  double* data_out = LOOP_MALLOC(double, ngen_elems * ncomps);
  for (unsigned i = 0; i < nedges; ++i) {
    if (gen_offset_of_edges[i] ==
        gen_offset_of_edges[i + 1])
      continue;
    double sum[MAX_NCOMPS] = {0};
    {
      unsigned f = elems_of_edges_offsets[i];
      unsigned e = elems_of_edges_offsets[i + 1];
      for (unsigned j = f; j < e; ++j) {
        unsigned elem = elems_of_edges[j];
        add_vectors(sum, data_in + elem * ncomps, sum, ncomps);
      }
    }
    {
      unsigned f = gen_offset_of_edges[i];
      unsigned e = gen_offset_of_edges[i + 1];
      double new_cavity_size = 0;
      for (unsigned j = f; j < e; ++j) {
        unsigned new_elem = nsame_elems + j;
        new_cavity_size += new_elem_sizes[new_elem];
      }
      for (unsigned j = f; j < e; ++j) {
        unsigned new_elem = nsame_elems + j;
        scale_vector(sum, new_elem_sizes[new_elem] / new_cavity_size,
            data_out + j * ncomps, ncomps);
      }
    }
  }
  return data_out;
}

static void swap_conserve_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_edges,
    unsigned const* offset_of_same_elems,
    double const* new_elem_sizes,
    struct const_tag* t)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned nedges = mesh_count(m, 1);
  double* same_data = doubles_expand(nelems, t->ncomps,
      t->d.f64, offset_of_same_elems);
  unsigned nsame_elems = uints_at(offset_of_same_elems, nelems);
  double* gen_data = swap_conserve_data(m, gen_offset_of_edges,
      nsame_elems, new_elem_sizes, t);
  double* data_out = concat_doubles(t->ncomps,
      same_data, nsame_elems,
      gen_data, uints_at(gen_offset_of_edges, nedges));
  loop_free(same_data);
  loop_free(gen_data);
  mesh_add_tag(m_out, elem_dim, TAG_F64, t->name, t->ncomps, data_out);
}

/* TODO: try to consolidate with coarsen_conserve */

void swap_conserve(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_edges,
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
    if (t->type == TAG_F64)
      swap_conserve_tag(m, m_out, gen_offset_of_edges,
          offset_of_same_elems, new_elem_sizes, t);
  }
  loop_free(new_elem_sizes);
}
