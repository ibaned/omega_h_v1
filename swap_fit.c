#include "swap_fit.h"

#include <assert.h>

#include "algebra.h"
#include "arrays.h"
#include "element_field.h"
#include "loop.h"
#include "mesh.h"
#include "qr.h"
#include "size.h"

#define MAX_NCOMPS 32

LOOP_KERNEL(swap_fit_cavity,
    unsigned ncomps,
    unsigned nsame_elems,
    unsigned const* gen_offset_of_edges,
    unsigned const* elems_of_edges_offsets,
    unsigned const* elems_of_edges,
    double const* old_elem_coords,
    double const* new_elem_coords,
    double const* data_in,
    double* gen_data)

  if (gen_offset_of_edges[i] ==
      gen_offset_of_edges[i + 1])
    return;
  unsigned f = elems_of_edges_offsets[i];
  unsigned e = elems_of_edges_offsets[i + 1];
  unsigned npts = e - f;
  if (npts > MAX_PTS)
    npts = MAX_PTS;
  enum { FIT, AVERAGE } method = FIT;
  if (npts < 4)
    method = AVERAGE;
  double q[MAX_PTS][MAX_PTS];
  double r[MAX_PTS][4];
  if (method == FIT) {
    double a[MAX_PTS][4];
    for (unsigned j = 0; j < npts; ++j) {
      unsigned old_elem = elems_of_edges[f + j];
      a[j][0] = 1.0;
      copy_vector(old_elem_coords + 3 * old_elem,
          &(a[j][1]), 3);
    }
    unsigned rank = qr_decomp2(a, q, r, npts);
    if (rank < 4)
      method = AVERAGE;
  }
  for (unsigned j = 0; j < ncomps; ++j) {
    double b[MAX_PTS];
    double c[4];
    double avg = 0;
    if (method == FIT) {
      for (unsigned k = 0; k < npts; ++k) {
        unsigned old_elem = elems_of_edges[f + k];
        b[k] = data_in[old_elem * ncomps + j];
      }
      qr_solve2(q, r, b, npts, c);
    } else {
      avg = 0;
      for (unsigned k = f; k < e; ++k) {
        unsigned old_elem = elems_of_edges[k];
        avg += data_in[old_elem * ncomps + j];
      }
      avg /= (e - f);
    }
    unsigned nf = gen_offset_of_edges[i];
    unsigned ne = gen_offset_of_edges[i + 1];
    for (unsigned gen_elem = nf; gen_elem < ne; ++gen_elem) {
      unsigned new_elem = gen_elem + nsame_elems;
      if (method == FIT)
        gen_data[gen_elem * ncomps + j] =
          c[0] + dot_product(c + 1, new_elem_coords + new_elem * 3, 3);
      else
        gen_data[gen_elem * ncomps + j] = avg;
    }
  }
}

static double* swap_fit_data(
    struct mesh* m,
    unsigned const* gen_offset_of_edges,
    unsigned nsame_elems,
    double const* old_elem_coords,
    double const* new_elem_coords,
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
  double* gen_data = LOOP_MALLOC(double, ngen_elems * ncomps);
  LOOP_EXEC(swap_fit_cavity, nedges,
      ncomps,
      nsame_elems,
      gen_offset_of_edges,
      elems_of_edges_offsets,
      elems_of_edges,
      old_elem_coords,
      new_elem_coords,
      data_in,
      gen_data);
  return gen_data;
}

static void swap_fit_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned const* gen_offset_of_edges,
    unsigned const* offset_of_same_elems,
    double const* old_elem_coords,
    double const* new_elem_coords,
    struct const_tag* t)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned nedges = mesh_count(m, 1);
  double* same_data = doubles_expand(nelems, t->ncomps,
      t->d.f64, offset_of_same_elems);
  unsigned nsame_elems = uints_at(offset_of_same_elems, nelems);
  double* gen_data = swap_fit_data(m, gen_offset_of_edges,
      nsame_elems, old_elem_coords, new_elem_coords, t);
  double* data_out = concat_doubles(t->ncomps,
      same_data, nsame_elems,
      gen_data, uints_at(gen_offset_of_edges, nedges));
  loop_free(same_data);
  loop_free(gen_data);
  add_tag2(mesh_tags(m_out, elem_dim), TAG_F64, t->name, t->ncomps,
      t->transfer_type, data_out);
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
    if (mesh_get_tag(m, elem_dim, i)->type == TAG_F64)
      break;
  if (i == mesh_count_tags(m, elem_dim))
    return;
  assert(elem_dim == 3);
  mesh_interp_to_elems(m, "coordinates");
  mesh_interp_to_elems(m_out, "coordinates");
  double const* old_elem_coords =
    mesh_find_tag(m, elem_dim, "coordinates")->d.f64;
  double const* new_elem_coords =
    mesh_find_tag(m_out, elem_dim, "coordinates")->d.f64;
  for (i = 0; i < mesh_count_tags(m, elem_dim); ++i) {
    struct const_tag* t = mesh_get_tag(m, elem_dim, i);
    if ((t->type == TAG_F64) && (t->transfer_type == OSH_TRANSFER_POINTWISE))
      swap_fit_tag(m, m_out, gen_offset_of_edges, offset_of_same_elems,
          old_elem_coords, new_elem_coords, t);
  }
}

