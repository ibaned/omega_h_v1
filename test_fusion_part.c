#include "algebra.h"
#include "bcast.h"
#include "comm.h"
#include "element_field.h"
#include "mesh.h"
#include "parallel_inertial_bisect.h"
#include "parallel_mesh.h"
#include "vtk.h"

int main()
{
  comm_init();
  struct mesh* m = 0;
  if (comm_rank() == 0)
    m = read_vtu("xgc.vtu");
  m = bcast_mesh_metadata(m);
  mesh_number_simply(m, 0);
  mesh_interp_to_elems(m, "coordinates");
  double* elem_coords = mesh_find_tag(m, 2, "coordinates")->d.f64;
  unsigned nelems = mesh_count(m, 2);
  double const center[2] = {1.75, 0};
  for (unsigned i = 0; i < nelems; ++i) {
    double rel[2];
    subtract_vectors(elem_coords + i * 3, center, rel, 2);
    elem_coords[i * 3 + 0] = vector_norm(rel, 2);
    elem_coords[i * 3 + 1] = 0;
    elem_coords[i * 3 + 2] = 0;
  }
  balance_mesh_inertial(&m);
  write_parallel_vtu(m, "split.pvtu");
  free_mesh(m);
  comm_fini();
}

