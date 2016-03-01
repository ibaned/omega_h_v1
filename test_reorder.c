#include "comm.h"
#include "mesh.h"
#include "reorder.h"
#include "vtk_io.h"

int main(int argc, char** argv)
{
  comm_init();
  (void) argc;
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_add_tag(m, 0, TAG_U32, "order", 1, compute_ordering(m));
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  comm_fini();
}
