#include "inertia.h"

#include "element_field.h"
#include "loop.h"
#include "mesh.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  (void) argc;
  struct mesh* m = read_vtk(argv[1]);
  mesh_interp_to_elems(m, "coordinates");
  unsigned n = mesh_count(m, mesh_dim(m));
  double const* coords = mesh_find_tag(m, mesh_dim(m), "coordinates")->data;
  double* masses = loop_malloc(sizeof(double) * n);
  for (unsigned i = 0; i < n; ++i)
    masses[i] = 1;
  unsigned* in;
  local_inertial_mark(n, coords, masses, 0.01, &in);
  mesh_add_tag(m, mesh_dim(m), TAG_U32, "part", 1, in);
  write_vtk(m, argv[2]);
  free_mesh(m);
}
