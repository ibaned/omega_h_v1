#include <cassert>
#include <cstdio>

#include "arrays.hpp"
#include "algebra.hpp"
#include "include/omega_h.h"
#include "mark.hpp"
#include "mesh.hpp"
#include "coarsen.hpp"
#include "vtk_io.hpp"
#include "adapt.hpp"
#include "loop.hpp"
#include "smooth.hpp"

using namespace omega_h;

#define MOTION 0

LOOP_KERNEL(move_object_vert,
    unsigned const* object_verts,
    double* coords)
  if (object_verts[i]) {
#if MOTION == 0
    coords[i * 3 + 0] += 0.02;
#elif MOTION == 1
    double x[3];
    double const mid[3] = {.5, .5, 0};
    subtract_vectors(coords + i * 3, mid, x, 3);
    double polar_a = atan2(x[1], x[0]);
    double polar_r = vector_norm(x, 3);
    double rot_a = PI / 16;
    double dest_a = polar_a + rot_a;
    x[0] = cos(dest_a) * polar_r;
    x[1] = sin(dest_a) * polar_r;
    x[2] = 0;
    add_vectors(x, mid, coords + i * 3, 3);
#elif MOTION == 2
    coords[i * 3 + 1] += 0.02;
#endif
  }
}

#if MOTION == 2
LOOP_KERNEL(move_object_vert2,
    unsigned const* object_verts,
    double* coords)
  if (object_verts[i]) {
    coords[i * 3 + 1] -= 0.02;
  }
}
#endif

static void move_mesh(struct mesh* m)
{
  unsigned nverts = mesh_count(m, 0);
  double* coords = mesh_find_tag(m, 0, "coordinates")->d.f64;
  unsigned* object_verts = mesh_mark_class_closure_verts(m, 3, 72);
  LOOP_EXEC(move_object_vert, nverts, object_verts, coords);
  loop_free(object_verts);
#if MOTION == 2
  unsigned* object_verts2 = mesh_mark_class_closure_verts(m, 2, 32);
  LOOP_EXEC(move_object_vert2, nverts, object_verts2, coords);
  loop_free(object_verts2);
#endif
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
  for (unsigned i = 0; i < 12; ++i) {
    move_mesh(m);
  //mesh_adapt(m, 0.3, 0.3, 4, 50);
    write_vtk_step(m);
  }
  free_mesh(m);
  osh_fini();
}

