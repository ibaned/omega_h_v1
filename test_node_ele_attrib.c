#include <assert.h>

#include "arrays.h"
#include "comm.h"
#include "node_ele_io.h"
#include "mesh.h"

int main(int argc, char** argv)
{
  assert(argc == 5);
  comm_init();
  struct mesh* m = read_dot_node(argv[1]);
  mesh_add_tag(m, 0, TAG_F64, "attributes", 3,
      doubles_filled(mesh_count(m, 0), 4.2));
  read_dot_ele(m, argv[2]);
  write_dot_node(m, argv[3]);
  write_dot_ele(m, argv[4]);
  free_mesh(m);
  comm_fini();
}

