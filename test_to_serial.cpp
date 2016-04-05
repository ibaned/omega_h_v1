#include <cassert>

#include "include/omega_h.h"
#include "mesh.hpp"
#include "comm.hpp"
#include "parallel_mesh.hpp"
#include "parallel_to_serial.hpp"
#include "vtk_io.hpp"

using namespace omega_h;

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 3);
  struct mesh* m = read_mesh_vtk(argv[1]);
  parallel_to_serial(m);
  if (!comm_rank())
    write_mesh_vtk(m, argv[2]);
  free_mesh(m);
  osh_fini();
}
