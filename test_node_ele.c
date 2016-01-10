#include <assert.h>

#include "comm.h"
#include "node_ele_io.h"
#include "mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 4);
  comm_init();
  struct mesh* m = read_dot_node(argv[1]);
  read_dot_ele(m, argv[2]);
  write_mesh_vtk(m, argv[3]);
  free_mesh(m);
  comm_fini();
}
