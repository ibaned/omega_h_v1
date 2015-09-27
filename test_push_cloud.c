#include "mesh.h"
#include "cloud.h"
#include "form_cloud.h"
#include "push_cloud.h"
#include "element_field.h"
#include "eval_field.h"
#include "vtk.h"
#include "algebra.h"

static void density_fun(double const* x, double* d)
{
  (void) x;
  *d = 2.1;
}

static void new_coord_fun(double const* x, double* nx)
{
  double const v[3] = {-.4, .4, 0};
  add_vectors(x, v, nx, 3);
}

int main()
{
  struct mesh* m = new_box_mesh(2);
  mesh_interp_to_elems(m, "coordinates");
  mesh_eval_field(m, 2, "cloud_density", 1, density_fun);
  struct cloud * c = form_cloud(m);
  write_vtk(m, "mesh.vtu");
  write_vtk_cloud(c, "cloud_0.vtu");
  cloud_eval_field(c, "new_coords", 3, new_coord_fun);
  write_vtk_cloud(c, "cloud_1.vtu");
  cloud_free_tag(c, "coordinates");
  cloud_rename_tag(c, "new_coords", "coordinates");
  write_vtk_cloud(c, "cloud_2.vtu");
  free_cloud(c);
  free_mesh(m);
}
