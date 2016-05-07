#include <cassert>

#include "arrays.hpp"
#include "include/omega_h.h"
#include "mark.hpp"
#include "mesh.hpp"
#include "coarsen.hpp"
#include "vtk_io.hpp"
#include "adapt.hpp"
#include "loop.hpp"
#include "smooth.hpp"

using namespace omega_h;

LOOP_KERNEL(move_object_vert,
    unsigned const* object_verts,
    double* coords)
  if (object_verts[i])
    coords[i * 3 + 0] += 0.02;
}

static void move_mesh(struct mesh* m)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned* object_verts = mesh_mark_class_closure_verts(m, 2, 21);
  double* coords = mesh_find_tag(m, 0, "coordinates")->d.f64;
  LOOP_EXEC(move_object_vert, nverts, object_verts, coords);
  loop_free(object_verts);
  mesh_smooth_field(m, "coordinates", 1e-4, 50);
}

int main(int argc, char** argv)
{
  assert(argc == 2);
  osh_init(&argc, &argv);
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_add_tag(m, 0, TAG_F64, "adapt_size", 1,
      filled_array(mesh_count(m, 0), 0.1));
  start_vtk_steps("imr");
  write_vtk_step(m);
  for (unsigned i = 0; i < 11; ++i) {
    move_mesh(m);
    mesh_adapt(m, 0.3, 0.3, 4, 50);
    write_vtk_step(m);
  }
  free_mesh(m);
  osh_fini();
}


