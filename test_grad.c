#include "mesh.h"
#include "eval_field.h"
#include "element_gradients.h"
#include "vtk.h"

static void fld_fun(double const x[3], double v[])
{
  v[0] = x[0] * 5 + x[1] * 2 + x[2] * -3;
}

int main()
{
  struct mesh* m = new_box_mesh(3);
  mesh_eval_field(m, "foo", 1, fld_fun);
  mesh_element_gradients(m, "foo");
  write_vtk(m, "grad.vtu");
  free_mesh(m);
}
