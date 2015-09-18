#include "vtk.h"
#include "mesh.h"

int main()
{
//{
//  struct mesh* m = new_box_mesh(3);
//  write_vtk(m, "one.vtu");
//  free_mesh(m);
//}
  {
    struct mesh* m = read_vtk("one.vtu");
    write_vtk(m, "two.vtu");
    free_mesh(m);
  }
}
