#include "node_ele_io.h"
#include "vtk.h"
#include "mesh.h"
#include <assert.h>

int main(int argc, char** argv)
{
  assert(argc == 4);
  struct mesh* m = read_dot_node(argv[1]);
  read_dot_ele(m, argv[2]);
  write_vtk(m, argv[3]);
}
