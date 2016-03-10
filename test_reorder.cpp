#include "include/omega_h.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "reorder.hpp"
#include "shuffle_mesh.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  (void) argc;
  struct mesh* m = read_mesh_vtk(argv[1]);
  unsigned* vert_num = compute_ordering(m);
  shuffle_mesh(m, vert_num);
  loop_free(vert_num);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}
