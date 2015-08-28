#include "mesh.h"
#include "eval_field.h"
#include "element_gradients.h"
#include "refine_by_size.h"
#include "recover_by_volume.h"
#include "algebra.h"
#include "vtk.h"
#include "size_from_hessian.h"

static void dye_fun(double const coords[3], double v[])
{
  double x[3];
  double const l[3] = {.25, .5, 0};
  double const r[3] = {.75, .5, 0};
  double dir = 1;
  subtract_vectors(coords, l, x, 3);
  if (vector_norm(x, 3) > .25) {
    dir = -1;
    subtract_vectors(coords, r, x, 3);
  }
  if (vector_norm(x, 3) > .25) {
    v[0] = 0;
    return;
  }
  v[0] = 4 * dir * (.25 - vector_norm(x, 3));
}

static void size_fun(double const x[], double s[])
{
  s[0] = 0.1;
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  mesh_eval_field(m, "adapt_size", 1, size_fun);
  while (refine_by_size(&m));
  mesh_eval_field(m, "dye", 1, dye_fun);
  mesh_element_gradients(m, "dye");
  mesh_recover_by_volume(m, "grad_dye");
  mesh_element_gradients(m, "rcov_grad_dye");
  mesh_recover_by_volume(m, "grad_rcov_grad_dye");
  mesh_free_nodal_field(m, "adapt_size");
  double weight = 0.1 / 75.0;
  mesh_size_from_hessian(m, "rcov_grad_rcov_grad_dye", &weight, 0.1, 0.2);
  write_vtk(m, "grad.vtu");
  free_mesh(m);
}
