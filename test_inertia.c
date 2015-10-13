#include "element_field.h"
#include "inertia.h"
#include "mesh.h"
#include "tag.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  (void) argc;
  struct mesh* m = read_vtu(argv[1]);
  mesh_interp_to_elems(m, "coordinates");
  unsigned n = mesh_count(m, mesh_dim(m));
  double const* coords = mesh_find_tag(m, mesh_dim(m), "coordinates")->d.f64;
  double* masses = 0;
  unsigned* in = local_inertia_mark(n, coords, masses);
  mesh_add_tag(m, mesh_dim(m), TAG_U32, "part", 1, in);
  write_vtu(m, argv[2]);
  free_mesh(m);
}
