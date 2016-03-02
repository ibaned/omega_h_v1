#include "comm.h"
#include "mesh.h"
#include "reorder.h"
#include "shuffle_mesh.h"
#include "vtk_io.h"

int main(int argc, char** argv)
{
  comm_init();
  (void) argc;
  struct mesh* m = read_mesh_vtk(argv[1]);
  unsigned* vert_num = compute_ordering(m);
  mesh_add_tag(m, 0, TAG_U32, "order", 1, vert_num);
  mesh_add_tag(m, 3, TAG_U32, "order", 1, number_ents(m, 3, vert_num));
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  comm_fini();
}
