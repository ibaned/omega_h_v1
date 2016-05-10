#include <cassert>
#include <cstdio>

#include "include/omega_h.h"
#include "mesh.hpp"
#include "refine.hpp"
#include "vtk_io.hpp"
#include "timer.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  assert(argc == 3);
  osh_init(&argc, &argv);
  struct mesh* m = read_mesh_vtk(argv[1]);
  Now t0 = now();
  uniformly_refine(m);
  Now t1 = now();
  printf("uniform refinement took %e seconds\n", seconds_between(t0,t1));
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}
