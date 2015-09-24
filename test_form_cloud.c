#include "mesh.h"
#include "form_cloud.h"
#include "vtk.h"
#include "loop.h"
#include "cloud.h"

int main()
{
  struct mesh* m = new_box_mesh(2);
  double* density = loop_malloc(sizeof(double) * 2);
  density[0] = 6.1;
  density[1] = 6.1;
  mesh_add_tag(m, 2, TAG_F64, "cloud_density", 1, density);
  struct cloud* c = form_cloud(m);
  write_vtk(m, "cloud_mesh.vtu");
  write_vtk_cloud(c, "cloud.vtu");
  free_cloud(c);
  free_mesh(m);
}
