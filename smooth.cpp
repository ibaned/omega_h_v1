#include "smooth.hpp"

#include "algebra.hpp"
#include "loop.hpp"
#include "mark.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"

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

static void smooth_field_iter(
    unsigned n,
    unsigned ncomps,
    unsigned const* interior,
    unsigned const* star_offsets,
    unsigned const* star,
    double** p_data)
{
  double* data_in = *p_data;
  double* data_out = LOOP_MALLOC(double, n * ncomps);
  LOOP_EXEC(smooth_field_vert, n, ncomps, interior, star_offsets, star,
      data_in, data_out);
  loop_free(data_in);
  *p_data = data_out;
}

unsigned mesh_smooth_field(struct mesh* m, char const* name,
    double tol, unsigned maxiter)
{
  if (mesh_is_paralle(m))
    mesh_ensure_ghosting(m, 1);
  unsigned n = mesh_count(m, 0);
  struct const_tag* t = mesh_find_tag(m, 0, name);
  unsigned ncomps = t->ncomps;
  double* data = doubles_copy(t->d.f64, n * ncomps);
  mesh_free_tag(m, 0, name);
  unsigned const* star_offsets = mesh_ask_star(m, 0, 1)->offsets;
  unsigned const* star = mesh_ask_star(m, 0, 1)->adj;
  unsigned* interior = mesh_mark_class(m, 0, mesh_dim(m), INVALID);
  unsigned i;
  for (i = 0; i < maxiter; ++i) {
    smooth_field_iter(n, ncomps, interior, star_offsets, star, &data);
    mesh_conform_doubles(m, 0, ncomps, &data);
  }
  loop_free(interior);
  mesh_add_tag(m, 0, TAG_F64, name, ncomps, data);
  return i;
}

}
