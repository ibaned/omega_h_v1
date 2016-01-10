#include "element_field.h"
#include "inertia.h"
#include "mesh.h"
#include "tag.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  (void) argc;
  struct mesh* m = read_mesh_vtk(argv[1]);
  mesh_interp_to_elems(m, "coordinates");
  unsigned n = mesh_count(m, mesh_dim(m));
  double const* coords = mesh_find_tag(m, mesh_dim(m), "coordinates")->d.f64;
  double* masses = 0;
  unsigned* in = mark_inertial_bisection(n, coords, masses, 0);
  mesh_add_tag(m, mesh_dim(m), TAG_U32, "part", 1, in);
  write_mesh_vtk(m, argv[2]);
  free_mesh(m);
}
