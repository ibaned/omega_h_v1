#include <assert.h>

#include "arrays.hpp"
#include "comm.hpp"
#include "node_ele_io.hpp"
#include "mesh.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  assert(argc == 4);
  comm_init();
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_add_tag(m, 0, TAG_F64, "attributes", 3,
      doubles_filled(mesh_count(m, 0) * 3, 4.2));
  write_dot_node(m, argv[2]);
  write_dot_ele(m, argv[3]);
  free_mesh(m);
  comm_fini();
}

