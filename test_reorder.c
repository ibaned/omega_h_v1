#include "comm.h"
#include "loop.h"
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
  shuffle_mesh(m, vert_num);
  loop_free(vert_num);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  comm_fini();
}
