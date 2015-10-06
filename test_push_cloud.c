#include "mesh.h"
#include "cloud.h"
#include "push_cloud.h"
#include "eval_field.h"
#include "vtk.h"
#include "algebra.h"

static double const* cen;

static void new_coord_fun(double const* x, double* nx)
{
  double v[3];
  subtract_vectors(x, cen, v, 3);
  double r = vector_norm(v, 3);
  double a = atan2(v[1], v[0]);
  a += PI / 8;
  v[0] = cos(a) * r;
  v[1] = sin(a) * r;
  v[2] = 0;
  add_vectors(v, cen, nx, 3);
}

int main()
{
  struct mesh* m = read_vtk("xgc.vtu");
  /* assume first vertex is center vertex */
  cen = mesh_find_tag(m, 0, "coordinates")->data;
  struct cloud* c = read_vtk_cloud("cloud.vtu");
  cloud_eval_field(c, "new_coords", 3, new_coord_fun);
  cloud_free_tag(c, "coordinates");
  cloud_rename_tag(c, "new_coords", "coordinates");
  push_cloud(c, m);
  write_vtk_cloud(c, "pushed_cloud.vtu");
  free_cloud(c);
  free_mesh(m);
}
