#include "mesh.h"
#include "vtk.h"

int main()
{
  struct mesh* m = new_box_mesh(3);
  write_vtu(m, "test.vtu");
  free_mesh(m);
}
