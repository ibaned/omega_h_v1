#include "mesh.h"
#include "form_cloud.h"
#include "vtk.h"
#include "loop.h"
#include "cloud.h"
#include "element_field.h"
#include "eval_field.h"
#include "algebra.h"

static void density_fun(double const x[3], double* d)
{
  static double c[3] = {0.9, 0, 0};
  double r = vector_distance(x, c, 3);
  if (0.3 < r && r < 0.45)
    *d = 1000;
  else
    *d = 0;
}

int main()
{
  struct mesh* m = read_vtk("xgc.vtu");
  mesh_interp_to_elems(m, "coordinates");
  mesh_eval_field(m, 2, "cloud_density", 1, density_fun);
  struct cloud* c = form_cloud(m);
  write_vtk(m, "cloud_mesh.vtu");
  write_vtk_cloud(c, "cloud.vtu");
  free_cloud(c);
  free_mesh(m);
}
