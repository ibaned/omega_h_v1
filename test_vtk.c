#include "vtk.h"
#include "tables.h"

int main()
{
  unsigned dim = 3;
  write_vtk("test.vtu",
      dim,
      the_box_nelems[dim],
      the_box_nverts[dim],
      the_box_conns[dim],
      the_box_coords[dim],
      0);
}
