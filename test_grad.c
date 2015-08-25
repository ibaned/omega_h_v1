#include "mesh.h"
#include "element_gradients.h"
#include "vtk.h"

int main()
{
  struct mesh* m = new_box_mesh(2);
  mesh_element_gradients(m, "coordinates");
  write_vtk(m, "grad.vtu");
}
