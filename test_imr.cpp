#include <cassert>

#include "arrays.hpp"
#include "include/omega_h.h"
#include "mesh.hpp"
#include "coarsen.hpp"
#include "vtk_io.hpp"
#include "adapt.hpp"
#include "loop.hpp"
#include "smooth.hpp"

using namespace omega_h;

LOOP_KERNEL(move_vert,
    unsigned const* class_id,
    double* coords)
  switch (class_id[i]) {
    case 12:
    case 13:
    case 14:
    case 15:
    case 17:
    case 18:
    case 19:
    case 20:
    case 21:
      coords[i * 3 + 0] += 0.02;
  }
}

static void move_mesh(struct mesh* m)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned const* class_id = mesh_find_tag(m, 0, "class_id")->d.u32;
  double* coords = mesh_find_tag(m, 0, "coordinates")->d.f64;
  LOOP_EXEC(move_vert, nverts, class_id, coords);
  mesh_smooth_field(m, "coordinates", 0.01, 50);
}

int main(int argc, char** argv)
{
  assert(argc == 2);
  osh_init(&argc, &argv);
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_add_tag(m, 0, TAG_F64, "adapt_size", 1,
      filled_array(mesh_count(m, 0), 0.1));
  start_vtk_steps("imr");
  for (unsigned i = 0; i < 1; ++i) {
    move_mesh(m);
  //mesh_adapt(m, 0.3, 0.3, 4, 50);
    write_vtk_step(m);
  }
  free_mesh(m);
  osh_fini();
}


