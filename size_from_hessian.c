#include "size_from_hessian.h"
#include "algebra.h"
#include "doubles.h"
#include <stdlib.h>
#include <assert.h>

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
  double* out = malloc(sizeof(double) * nverts);
  unsigned nsol_comps = nhess_comps / 9;
  for (unsigned i = 0; i < nverts; ++i) {
    double const* hess = hessians + i * nhess_comps;
    double norm_sq = 0;
    for (unsigned j = 0; j < nsol_comps; ++j) {
      double hess_norm_sq = dot_product(hess, hess, 9);
      double hess_w = 1;
      if (sol_comp_weights)
        hess_w = sol_comp_weights[j];
      norm_sq += hess_w + hess_norm_sq;
      hess += 9;
    }
    out[i] = sqrt(norm_sq);
  }
  double min_norm = doubles_min(out, nverts);
  double max_norm = doubles_max(out, nverts);
  double a = 0;
  if (max_norm != min_norm)
    a = (max_h - min_h) / (max_norm - min_norm);
  double b = min_h - min_norm;
  for (unsigned i = 0; i < nverts; ++i)
    out[i] = a * out[i] + b;
  return out;
}
