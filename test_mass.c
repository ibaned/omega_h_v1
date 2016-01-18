#include "adapt.h"
#include "comm.h"
#include "element_field.h"
#include "loop.h"
#include "mesh.h"
#include "size.h"
#include "vtk.h"

static double const good_qual_floor = 0.3;
static double const size_floor = 1. / 3.;
static unsigned const nsliver_layers = 4;
static unsigned const max_ops = 50;

int main()
{
  comm_init();
  struct mesh* m = read_mesh_vtk("before_10.vtu");
  unsigned dim = mesh_dim(m);
  if (mesh_find_tag(m, dim, "mass"))
    mesh_free_tag(m, dim, "mass");
  { //set mass field to test conservative transfer
    mesh_interp_to_elems(m, "coordinates");
    double const* ecoords = mesh_find_tag(m, dim, "coordinates")->d.f64;
    double* mass = mesh_element_sizes(m);
    for (unsigned i = 0; i < mesh_count(m, dim); ++i)
      if (ecoords[i * 3 + 0] < 0.5)
        mass[i] = 0;
    mesh_free_tag(m, mesh_dim(m), "coordinates");
    mesh_add_tag(m, dim, TAG_F64, "mass", 1, mass);
  }
  write_mesh_vtk(m, "output_0.vtu");
  mesh_adapt(&m, size_floor, good_qual_floor, nsliver_layers, max_ops);
  write_mesh_vtk(m, "output_1.vtu");
  free_mesh(m);
  comm_fini();
}
