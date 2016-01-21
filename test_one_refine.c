#include <assert.h>

#include "comm.h"
#include "mesh.h"
#include "refine_by_size.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 3);
  comm_init();
  struct mesh* m = read_mesh_vtk(argv[1]);
  uniformly_refine(&m);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  comm_fini();
}
