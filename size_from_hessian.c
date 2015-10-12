#include "size_from_hessian.h"

#include <assert.h>

#include "algebra.h"
#include "loop.h"
#include "mesh.h"
#include "tag.h"

double* size_from_hessian(
    unsigned nverts,
    unsigned nhess_comps,
    double const* hessians,
    double const* sol_comp_weights,
    double min_h,
    double max_h)
{
  assert(nhess_comps % 9 == 0);
  assert(max_h > min_h);
  assert(min_h > 0);
  double* out = loop_malloc(sizeof(double) * nverts);
  unsigned nsol_comps = nhess_comps / 9;
  for (unsigned i = 0; i < nverts; ++i) {
    double const* hess = hessians + i * nhess_comps;
    double total = 0;
    for (unsigned j = 0; j < nsol_comps; ++j) {
      double hess_norm = vector_norm(hess, 9);
      double hess_w = 1;
      if (sol_comp_weights)
        hess_w = sol_comp_weights[j];
      assert(hess_w >= 0);
      total += hess_w * hess_norm;
      hess += 9;
    }
    out[i] = total;
  }
  for (unsigned i = 0; i < nverts; ++i) {
    out[i] = max_h - out[i];
    if (out[i] < min_h)
      out[i] = min_h;
  }
  return out;
}

struct const_tag* mesh_size_from_hessian(struct mesh* m, char const* hess_name,
    double const* sol_comp_weights, double min_h, double max_h)
{
  struct const_tag* hf = mesh_find_tag(m, 0, hess_name);
  double* data = size_from_hessian(mesh_count(m, 0),
      hf->ncomps, hf->data, sol_comp_weights, min_h, max_h);
  return mesh_add_tag(m, 0, TAG_F64, "adapt_size", 1, data);
}
