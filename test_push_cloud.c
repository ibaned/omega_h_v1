#include "mesh.h"
#include "cloud.h"
#include "form_cloud.h"
#include "push_cloud.h"
#include "element_field.h"
#include "eval_field.h"
#include "vtk.h"

static void density_fun(double const x[3], double* d)
{
  (void) x;
  *d = 2.1;
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  mesh_interp_to_elems(m, "coordinates");
  mesh_eval_field(m, 2, "cloud_density", 1, density_fun);
  struct cloud * c = form_cloud(m);
  write_vtk(m, "mesh.vtu");
  write_vtk_cloud(c, "cloud_0.vtu");
  free_cloud(c);
  free_mesh(m);
}
