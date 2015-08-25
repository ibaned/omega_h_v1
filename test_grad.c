#include "mesh.h"
#include "eval_field.h"
#include "element_gradients.h"
#include "refine_by_size.h"
#include "vtk.h"

static void fld_fun(double const x[3], double v[])
{
  v[0] = x[0] * x[0];
}

static double size_fun(double const x[])
{
  return 0.1;
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  while (refine_by_size(&m, size_fun));
  mesh_eval_field(m, "foo", 1, fld_fun);
  mesh_element_gradients(m, "foo");
  write_vtk(m, "grad.vtu");
  free_mesh(m);
}
