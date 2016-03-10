#include <assert.h>

#include "include/omega_h.hpp"
#include "node_ele_io.hpp"
#include "mesh.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 4);
  struct mesh* m = read_dot_node(argv[1]);
  read_dot_ele(m, argv[2]);
  write_mesh_vtk(m, argv[3]);
  free_mesh(m);
  osh_fini();
}
