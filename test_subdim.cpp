#include <assert.h>
#include <stdlib.h>

#include "comm.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "subset.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  comm_init();
  assert(argc == 4);
  struct mesh* m = read_mesh_vtk(argv[1]);
  unsigned dim = (unsigned) atoi(argv[2]);
  assert(dim <= 3);
  if (mesh_is_parallel(m))
    for (unsigned d = 0; d <= mesh_dim(m); ++d)
      if (mesh_has_dim(m, d))
        mesh_parallel_to_tags(m, d);
  unsigned* offsets = uints_linear(mesh_count(m, dim) + 1, 1);
  struct mesh* sm = subset_mesh(m, dim, offsets);
  loop_free(offsets);
  free_mesh(m);
  write_mesh_vtk(sm, argv[3]);
  free_mesh(sm);
  comm_fini();
}
