#include <cassert>
#include <cstdio>

#include "include/omega_h.h"
#include "mesh.hpp"
#include "refine.hpp"
#include "vtk_io.hpp"
#include "timer.hpp"
#include "quality.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  assert(argc == 3);
  osh_init(&argc, &argv);
  struct mesh* m = read_mesh_vtk(argv[1]);
  printf("initial nelems %u, initial minqual %e\n",
      mesh_count(m, mesh_dim(m)), mesh_min_quality(m));
  Now t0 = now();
  uniformly_refine(m, 0.3);
  Now t1 = now();
  printf("%f seconds, %u elements, minqual %e\n", seconds_between(t0,t1),
      mesh_count(m, mesh_dim(m)),
      mesh_min_quality(m));
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}
