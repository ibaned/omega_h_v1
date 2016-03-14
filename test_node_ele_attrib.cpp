#include <cassert>

#include "arrays.hpp"
#include "include/omega_h.h"
#include "node_ele_io.hpp"
#include "mesh.hpp"
#include "vtk_io.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 4);
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_add_tag(m, 0, TAG_F64, "attributes", 3,
      filled_array(mesh_count(m, 0) * 3, 4.2));
  write_dot_node(m, argv[2]);
  write_dot_ele(m, argv[3]);
  free_mesh(m);
  osh_fini();
}

