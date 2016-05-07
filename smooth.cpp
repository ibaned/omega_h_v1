#include "smooth.hpp"

#include <cstdio>

#include "algebra.hpp"
#include "arrays.hpp"
#include "comm.hpp"
#include "ghost_mesh.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mark.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "tables.hpp"

namespace omega_h {

LOOP_KERNEL(smooth_field_vert,
    unsigned ncomps,
    unsigned const* interior,
    unsigned const* star_offsets,
    unsigned const* star,
    double const* data_in,
    double* data_out)
  if (!interior[i]) {
    copy_vector(data_in + i * ncomps, data_out + i * ncomps, ncomps);
    return;
  }
  zero_vector(data_out + i * ncomps, ncomps);
  unsigned a = star_offsets[i];
  unsigned b = star_offsets[i + 1];
  for (unsigned j = a; j < b; ++j) {
    unsigned other = star[j];
    add_vectors(data_in + other * ncomps,
        data_out + i * ncomps,
        data_out + i * ncomps,
        ncomps);
  }
  scale_vector(data_out + i * ncomps, 1.0 / (b - a),
      data_out + i * ncomps, ncomps);
}

LOOP_KERNEL(check_diff,
    double const* a,
    double const* b,
    double tol,
    unsigned* diffs)
  diffs[i] = (fabs(a[i] - b[i]) > tol);
}

static unsigned smooth_field_iter(
    unsigned n,
    unsigned ncomps,
    unsigned const* interior,
    unsigned const* star_offsets,
    unsigned const* star,
    double tol,
    double** p_data)
{
  double* data_in = *p_data;
  double* data_out = LOOP_MALLOC(double, n * ncomps);
  LOOP_EXEC(smooth_field_vert, n, ncomps, interior, star_offsets, star,
      data_in, data_out);
  unsigned* diffs = LOOP_MALLOC(unsigned, n * ncomps);
  LOOP_EXEC(check_diff, n * ncomps, data_in, data_out, tol, diffs);
  loop_free(data_in);
  *p_data = data_out;
  unsigned not_done = uints_max(diffs, n * ncomps);
  loop_free(diffs);
  not_done = comm_max_uint(not_done);
  return not_done;
}

unsigned mesh_smooth_field(struct mesh* m, char const* name,
    double tol, unsigned maxiter)
{
  if (mesh_is_parallel(m))
    mesh_ensure_ghosting(m, 1);
  unsigned n = mesh_count(m, 0);
  struct const_tag* t = mesh_find_tag(m, 0, name);
  unsigned ncomps = t->ncomps;
  double* data = copy_array(t->d.f64, n * ncomps);
  mesh_free_tag(m, 0, name);
  unsigned const* star_offsets = mesh_ask_star(m, 0, 1)->offsets;
  unsigned const* star = mesh_ask_star(m, 0, 1)->adj;
  unsigned* interior = mesh_mark_class(m, 0, mesh_dim(m), INVALID);
  unsigned i;
  for (i = 0; i < maxiter; ++i) {
    unsigned not_done = smooth_field_iter(
        n, ncomps, interior, star_offsets, star, tol, &data);
    mesh_conform_array(m, 0, ncomps, &data);
    if (!not_done)
      break;
  }
  loop_free(interior);
  mesh_add_tag(m, 0, TAG_F64, name, ncomps, data);
  if (i == maxiter)
    printf("smoothing failed to converge after %u iterations\n", i);
  else
    printf("smoothing converged in %u iterations\n", i);
  return i;
}

}
