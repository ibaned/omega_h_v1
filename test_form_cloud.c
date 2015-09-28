#include "mesh.h"
#include "form_cloud.h"
#include "vtk.h"
#include "loop.h"
#include "cloud.h"
#include "element_field.h"
#include "eval_field.h"
#include "algebra.h"

static double const* cen;

static void density_fun(double const x[3], double* d)
{
  double r = vector_distance(x, cen, 3);
  if (0.1 < r && r < 0.15)
    *d = 1000;
  else
    *d = 0;
}

int main()
{
  struct mesh* m = read_vtk("xgc.vtu");
  /* assume first vertex is center vertex */
  cen = mesh_find_tag(m, 0, "coordinates")->data;
  mesh_interp_to_elems(m, "coordinates");
  mesh_eval_field(m, 2, "cloud_density", 1, density_fun);
  struct cloud* c = form_cloud(m);
  write_vtk(m, "cloud_mesh.vtu");
  write_vtk_cloud(c, "cloud.vtu");
  free_cloud(c);
  free_mesh(m);
}
